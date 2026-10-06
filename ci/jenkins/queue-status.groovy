// Install reviewed code as $JENKINS_HOME/init.groovy.d/gkeyll-queue-status.groovy.
// This runs on the controller, never in a candidate checkout or build executor.
import com.cloudbees.plugins.credentials.CredentialsProvider
import com.cloudbees.plugins.credentials.common.StandardUsernamePasswordCredentials
import com.cloudbees.plugins.credentials.domains.URIRequirementBuilder
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import hudson.ExtensionList
import hudson.model.Job
import hudson.model.ParametersAction
import hudson.model.ParametersDefinitionProperty
import hudson.model.Queue
import hudson.model.Result
import hudson.model.Run
import hudson.model.StringParameterValue
import hudson.model.TaskListener
import hudson.model.listeners.RunListener
import hudson.model.queue.CauseOfBlockage
import hudson.model.queue.QueueListener
import hudson.model.queue.QueueTaskDispatcher
import hudson.slaves.EnvironmentVariablesNodeProperty
import hudson.security.ACL
import hudson.util.AtomicFileWriter
import jenkins.model.Jenkins
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.Executors
import java.util.concurrent.ThreadFactory
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicBoolean
import java.util.logging.Logger

class GkeyllQueueStatus {
    static final Logger LOG = Logger.getLogger('gkeyll.queue-status')
    static final List INTERNAL = ['CI_QUEUE_ID', 'CI_QUEUE_COMMIT', 'CI_QUEUE_CONTEXT',
        'CI_QUEUE_STARTED', 'CI_QUEUE_CANCEL_DESCRIPTION']
    static final Map PLATFORMS = [
        personal: ['PERSONAL', 'gkeyll-ci-personal', null],
        stellar_cpu: ['STELLAR_CPU', 'gkeyll-ci-stellar_cpu', 'continuous-integration/jenkins/stellar_cpu'],
        perlmutter_gpu: ['PERLMUTTER_GPU', 'gkeyll-ci-perlmutter_gpu', 'continuous-integration/jenkins/perlmutter_gpu'],
        team: ['TEAM_WORKSTATION', 'gkeyll-ci-team-workstation', 'continuous-integration/jenkins/team-workstation']
    ]
    final Map<Long, Map> states = new ConcurrentHashMap<>()
    final File directory = new File(Jenkins.get().rootDir, 'gkeyll-queue-status')
    final AtomicBoolean checking = new AtomicBoolean(false)
    final def workers = Executors.newFixedThreadPool(4, { Runnable task ->
        def thread = new Thread(task, 'gkeyll-queue-status')
        thread.daemon = true
        thread
    } as ThreadFactory)
    final def timer = Executors.newSingleThreadScheduledExecutor({ Runnable task ->
        def thread = new Thread(task, 'gkeyll-queue-supersession')
        thread.daemon = true
        thread
    } as ThreadFactory)

    void start() {
        timer.scheduleWithFixedDelay({ checkSuperseded() } as Runnable, 60, 60, TimeUnit.SECONDS)
    }

    Map configuration(Job job) {
        def vars = new HashMap(System.getenv())
        vars.putAll(Jenkins.get().globalNodeProperties.get(EnvironmentVariablesNodeProperty)?.envVars ?: [:])
        for (entry in PLATFORMS) {
            def (prefix, defaultJob, defaultContext) = entry.value
            def root = vars["GKEYLL_CI_QUEUE_${prefix}_JOB"] ?: defaultJob
            if (job.fullName != root && !(entry.key == 'team' && job.parent.fullName == root)) continue
            // The multibranch parent is not a build job. Only its direct children qualify.
            if (entry.key == 'team' && job.parent.fullName != root) return null
            def context = (entry.key in ['personal', 'team'] ? vars["${prefix}_STATUS_CONTEXT"] : null) ?: defaultContext
            return [platform: entry.key, credential: vars["${prefix}_GITHUB_CREDENTIAL_ID"]?.trim(),
                context: context?.trim()]
        }
        null
    }

    static Map parameters(def item) {
        def values = [:]
        def job = item instanceof Run ? item.parent : item.task
        job.getProperty(ParametersDefinitionProperty)?.parameterDefinitions?.each {
            if (it.defaultParameterValue != null) values[it.name] = it.defaultParameterValue.value
        }
        item.getAction(ParametersAction)?.allParameters?.each { values[it.name] = it.value }
        values
    }

    static String queuedDescription(def id) { "Gkeyll CI queued (Jenkins queue #${id})." }

    Map queuedParameters(long id, Job job) {
        def record = states.get(id)
        if (record?.ready?.get() && record.metadata) return record.metadata
        def file = new File(directory, "${id}.json")
        if (!file.isFile()) return [:]
        def saved = new JsonSlurper().parse(file, 'UTF-8')
        saved.job == job.fullName && saved.parameters?.CI_QUEUE_ID == id.toString() ? saved.parameters : [:]
    }

    void saveParameters(long id, Job job, Map values) {
        if (!directory.isDirectory() && !directory.mkdirs() && !directory.isDirectory()) {
            throw new IOException('Cannot create queue metadata directory')
        }
        def writer = new AtomicFileWriter(new File(directory, "${id}.json"), 'UTF-8')
        try {
            writer.write(JsonOutput.toJson([job: job.fullName, parameters: values]))
            writer.commit()
        } finally {
            writer.abort()
        }
    }

    void forget(long id) {
        states.remove(id)
        new File(directory, "${id}.json").delete()
    }

    static boolean validRef(String ref) {
        ref && ref ==~ /[A-Za-z0-9][A-Za-z0-9._\/-]*/ && !ref.contains('..') &&
            !ref.contains('//') && !ref.endsWith('/') && !ref.endsWith('.') && !ref.endsWith('.lock')
    }

    // Keep URLs fixed to GitHub and look up the credential in the job's scope.
    // Never log response bodies, request headers, or exception messages containing them.
    def request(Job job, Map config, String method, String path, Map body = null) {
        def credential = CredentialsProvider.lookupCredentialsInItem(
            StandardUsernamePasswordCredentials, job, ACL.SYSTEM2,
            URIRequirementBuilder.fromUri('https://api.github.com').build()
        ).find { it.id == config.credential }
        if (!credential) throw new IOException('GitHub credential is unavailable')
        def connection = new URL("https://api.github.com/repos/gkeyllorg/gkeyll/${path}").openConnection()
        connection.connectTimeout = 10000
        connection.readTimeout = 15000
        connection.instanceFollowRedirects = false
        connection.requestMethod = method
        connection.setRequestProperty('Accept', 'application/vnd.github+json')
        connection.setRequestProperty('Authorization', "Bearer ${credential.password.plainText}")
        connection.setRequestProperty('User-Agent', 'gkeyll-jenkins-queue-status')
        try {
            if (body != null) {
                connection.doOutput = true
                connection.setRequestProperty('Content-Type', 'application/json')
                connection.outputStream.withCloseable { it.write(JsonOutput.toJson(body).getBytes('UTF-8')) }
            }
            int status = connection.responseCode
            if (status < 200 || status >= 300) throw new IOException("GitHub HTTP ${status}")
            connection.inputStream.withCloseable { new JsonSlurper().parse(it, 'UTF-8') }
        } finally {
            connection.disconnect()
        }
    }

    Map selection(Job job, Map config, Map values) {
        def pr = values.CANDIDATE_PR?.toString()?.trim()
        def ref = values.CANDIDATE_REF?.toString()?.trim()
        def baseline = values.BASELINE_REF?.toString()?.trim()
        boolean automatic = config.platform == 'team' && !pr && !ref && !baseline
        if (automatic) {
            if (job.name == 'main') { ref = 'main'; baseline = 'main' }
            else if (job.name ==~ /PR-[1-9][0-9]*/) pr = job.name.substring(3)
            else return null
        }
        if (!!pr == !!ref || (pr && (baseline || !(pr ==~ /[1-9][0-9]*/)))) return null
        if (ref && (!validRef(ref) || !validRef(baseline))) return null
        String sha
        if (pr) {
            def pull = request(job, config, 'GET', "pulls/${pr}")
            if (automatic && pull.base.ref != 'main') return null
            sha = pull.head.sha
        } else if (ref ==~ /[0-9a-fA-F]{40}/) {
            sha = ref.toLowerCase()
        } else {
            sha = request(job, config, 'GET', "commits/${URLEncoder.encode(ref, 'UTF-8')}").sha
        }
        if (!(sha ==~ /[0-9a-f]{40}/)) throw new IOException('Invalid candidate SHA')
        [sha: sha, pr: pr ?: '', ref: ref ?: '', baseline: baseline ?: '']
    }

    Map begin(Queue.Item item) {
        if (!(item.task instanceof Job)) return null // Ignore Pipeline node() queue entries.
        def config = configuration(item.task)
        if (!config) return null
        def record = [ready: new AtomicBoolean(false), cancelled: new AtomicBoolean(false), config: config]
        def existing = states.putIfAbsent(item.id, record)
        if (existing != null) return existing
        workers.submit({
            synchronized (record) {
                try {
                    prepare(item, record)
                } catch (Exception failure) {
                    LOG.warning("Queue #${item.id}: GitHub queue reporting failed (${failure.class.simpleName}); build may proceed. Check credential, context and GitHub connectivity.")
                } finally {
                    record.ready.set(true)
                    Jenkins.get().queue.scheduleMaintenance()
                }
            }
            checkSuperseded()
        } as Runnable)
        record
    }

    void prepare(Queue.Item item, Map record) {
        def config = record.config
        if (record.cancelled.get()) return
        if (!config.credential || !config.context) throw new IOException('Configure GitHub credential and status context')
        def saved = queuedParameters(item.id, item.task)
        // Preserve the accepted SHA across a controller restart.
        def selected = saved.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/
            ? [sha: saved.CI_QUEUE_COMMIT, pr: saved.CANDIDATE_PR ?: '',
               ref: saved.CANDIDATE_REF ?: '', baseline: saved.BASELINE_REF ?: '']
            : selection(item.task, config, parameters(item))
        if (!selected || record.cancelled.get()) return
        def metadata = [CI_QUEUE_ID: item.id.toString(), CI_QUEUE_COMMIT: selected.sha,
            CI_QUEUE_CONTEXT: saved.CI_QUEUE_CONTEXT ?: config.context,
            CI_QUEUE_STARTED: 'false', CI_QUEUE_CANCEL_DESCRIPTION: '', CANDIDATE_PR: selected.pr,
            CANDIDATE_REF: selected.ref, BASELINE_REF: selected.baseline]
        // Do not change queue ParametersAction: doing so breaks Jenkins' duplicate
        // submission folding (e.g. every multibranch scan could enqueue another run).
        saveParameters(item.id, item.task, metadata)
        record.metadata = metadata
        if (record.cancelled.get()) return
        request(item.task, config, 'POST', "statuses/${selected.sha}",
            [state: 'pending', context: metadata.CI_QUEUE_CONTEXT, description: queuedDescription(item.id)])
    }

    void left(Queue.LeftItem item) {
        if (!item.cancelled) return
        def record = states.get(item.id)
        if (record != null) record.cancelled.set(true)
        workers.submit({
            // Wait for any in-flight pending POST before sending cancellation.
            synchronized (record ?: new Object()) {
                def metadata = record?.metadata ?: queuedParameters(item.id, item.task)
                finish(item.task, metadata, item.id, 'error',
                    record?.cancelDescription ?: 'Gkeyll CI cancelled while queued.')
                forget(item.id)
            }
        } as Runnable)
    }

    void started(Run run) {
        def metadata = queuedParameters(run.queueId, run.parent)
        // Clear internal parameters copied by a rebuild or supplied by a caller.
        def overrides = INTERNAL.collectEntries { [(it): ''] } + metadata
        def action = run.getAction(ParametersAction) ?: new ParametersAction([])
        def extra = new ParametersAction(overrides.collect { k, v -> new StringParameterValue(k, v) },
            INTERNAL + ['CANDIDATE_PR', 'CANDIDATE_REF', 'BASELINE_REF'])
        run.addOrReplaceAction(action.merge(extra))
        run.save()
        forget(run.queueId)
    }

    // A Pipeline itself starts on a flyweight executor, then its first node()
    // may wait in the queue. Include that wait, but never abort allocated work.
    static Run nodeRun(Queue.Item item) {
        def owner = item.task.ownerExecutable
        owner instanceof Run ? owner : null
    }

    void nodeLeft(Queue.LeftItem item, Run run) {
        def metadata = parameters(run)
        if (metadata.CI_QUEUE_ID != run.queueId.toString() || metadata.CI_QUEUE_STARTED == 'true') return
        if (!item.cancelled) {
            def update = new ParametersAction([new StringParameterValue('CI_QUEUE_STARTED', 'true')], INTERNAL)
            run.addOrReplaceAction(run.getAction(ParametersAction).merge(update))
            workers.submit({ run.save() } as Runnable)
        }
        // A cancelled node() completes as ABORTED; completed() closes its status.
    }

    void checkSuperseded() {
        if (!checking.compareAndSet(false, true)) return
        try {
            def candidates = []
            Jenkins.get().queue.items.each { item ->
                def run = item.task instanceof Job ? null : nodeRun(item)
                def job = run?.parent ?: (item.task instanceof Job ? item.task : null)
                def config = job ? configuration(job) : null
                if (!config) return
                def values = run ? parameters(run) : queuedParameters(item.id, job)
                if (!(values.CANDIDATE_PR ==~ /[1-9][0-9]*/) ||
                    !(values.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/) || values.CI_QUEUE_STARTED == 'true') return
                if (values.CI_QUEUE_ID != (run ? run.queueId : item.id).toString()) return
                candidates << [id: item.id, run: run, job: job, config: config, values: values]
            }
            // One API lookup per PR per sweep, even with many queued revisions.
            candidates.groupBy { it.values.CANDIDATE_PR }.each { pr, builds ->
                try {
                    def sample = builds.first()
                    def head = request(sample.job, sample.config, 'GET', "pulls/${pr}").head.sha
                    if (!(head ==~ /[0-9a-f]{40}/)) return
                    builds.findAll { it.values.CI_QUEUE_COMMIT != head }.each { candidate ->
                        cancelSuperseded(candidate, head)
                    }
                } catch (Exception failure) {
                    LOG.warning("PR #${pr}: cannot check queued revisions (${failure.class.simpleName}); leaving builds queued.")
                }
            }
        } catch (Exception failure) {
            LOG.warning("Cannot inspect the Jenkins queue (${failure.class.simpleName}); will retry next sweep.")
        } finally {
            checking.set(false)
        }
    }

    void cancelSuperseded(Map candidate, String head) {
        String description = "Gkeyll CI cancelled while queued; superseded by ${head.take(12)}."
        Queue.withLock({
            def live = Jenkins.get().queue.getItem(candidate.id as long)
            if (live == null || live instanceof Queue.LeftItem) return
            // Recheck under the scheduling lock: never cancel a node that has started.
            def current = candidate.run ? parameters(candidate.run) : queuedParameters(live.id, live.task)
            if (current.CI_QUEUE_COMMIT != candidate.values.CI_QUEUE_COMMIT || current.CI_QUEUE_STARTED == 'true') return
            if (candidate.run) {
                def update = new ParametersAction([
                    new StringParameterValue('CI_QUEUE_CANCEL_DESCRIPTION', description)], INTERNAL)
                candidate.run.addOrReplaceAction(candidate.run.getAction(ParametersAction).merge(update))
            } else {
                def record = states.get(candidate.id)
                if (record != null) record.cancelDescription = description
            }
            if (Jenkins.get().queue.cancel(live)) LOG.info("Queue #${candidate.id}: ${description}")
        } as Runnable)
    }

    void completed(Run run) {
        def metadata = parameters(run)
        if (metadata.CI_QUEUE_ID != run.queueId.toString()) return
        workers.submit({
            def state = run.result == Result.SUCCESS ? 'success' :
                (run.result in [Result.FAILURE, Result.UNSTABLE] ? 'failure' : 'error')
            def description = metadata.CI_QUEUE_CANCEL_DESCRIPTION ?:
                (run.result == Result.ABORTED && metadata.CI_QUEUE_STARTED != 'true'
                    ? 'Gkeyll CI cancelled while queued.' : "Gkeyll CI finished: ${run.result}.")
            finish(run.parent, metadata, run.queueId, state, description)
        } as Runnable)
    }

    void finish(Job job, Map metadata, long id, String state, String description) {
        if (metadata.CI_QUEUE_ID != id.toString() || !(metadata.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/)) return
        def config = configuration(job)
        if (!config?.credential) return
        try {
            def current = request(job, config, 'GET', "commits/${metadata.CI_QUEUE_COMMIT}/status?per_page=100")
                .statuses.find { it.context == metadata.CI_QUEUE_CONTEXT }
            // Preserve the Pipeline's detailed terminal status and another run's newer status.
            if (current?.state != 'pending' || current.description != queuedDescription(id)) return
            request(job, config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}",
                [state: state, context: metadata.CI_QUEUE_CONTEXT, description: description])
        } catch (Exception failure) {
            LOG.warning("Queue #${id}: GitHub final queue status failed (${failure.class.simpleName}).")
        }
    }
}

class GkeyllQueueListener extends QueueListener {
    final GkeyllQueueStatus service
    GkeyllQueueListener(GkeyllQueueStatus service) { this.service = service }
    @Override void onEnterWaiting(Queue.WaitingItem item) { service.begin(item) }
    @Override void onLeft(Queue.LeftItem item) {
        if (item.task instanceof Job) {
            if (service.configuration(item.task)) service.left(item)
        } else {
            def run = GkeyllQueueStatus.nodeRun(item)
            if (run && service.configuration(run.parent)) service.nodeLeft(item, run)
        }
    }
}

class GkeyllQueueDispatcher extends QueueTaskDispatcher {
    final GkeyllQueueStatus service
    GkeyllQueueDispatcher(GkeyllQueueStatus service) { this.service = service }
    @Override CauseOfBlockage canRun(Queue.Item item) {
        def record = service.begin(item)
        record != null && !record.ready.get() ? new CauseOfBlockage() {
            @Override String getShortDescription() { 'Reporting queued Gkeyll CI to GitHub' }
        } : null
    }
}

class GkeyllQueueRunListener extends RunListener<Run> {
    final GkeyllQueueStatus service
    GkeyllQueueRunListener(GkeyllQueueStatus service) { super(Run); this.service = service }
    @Override void onStarted(Run run, TaskListener listener) {
        if (service.configuration(run.parent)) service.started(run)
    }
    @Override void onCompleted(Run run, TaskListener listener) {
        if (service.configuration(run.parent)) service.completed(run)
    }
}

// Queue metadata is JSON; build metadata uses Jenkins' own ParametersAction.
// Neither requires custom classes during deserialization before this hook loads.
def service = new GkeyllQueueStatus()
ExtensionList.lookup(QueueListener).add(new GkeyllQueueListener(service))
ExtensionList.lookup(QueueTaskDispatcher).add(new GkeyllQueueDispatcher(service))
ExtensionList.lookup(RunListener).add(new GkeyllQueueRunListener(service))
Jenkins.get().queue.items.each { service.begin(it) }
service.start()
GkeyllQueueStatus.LOG.info('Gkeyll GitHub queue status listener installed')
return service
