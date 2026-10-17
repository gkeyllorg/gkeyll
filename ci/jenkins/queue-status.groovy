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
import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.time.Duration

class GkeyllQueueStatus {
    static final Logger LOG = Logger.getLogger('gkeyll.queue-status')
    static final List INTERNAL = ['CI_QUEUE_ID', 'CI_QUEUE_COMMIT', 'CI_QUEUE_CONTEXT',
        'CI_QUEUE_STARTED', 'CI_QUEUE_CANCEL_DESCRIPTION', 'CI_QUEUE_ENQUEUED_MS', 'CI_AGENT_STARTED_MS']
    static final HttpClient HTTP = HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(10)).build()
    static final Map PLATFORMS = [
        personal: ['PERSONAL', 'gkeyll-ci-personal', null],
        stellar_cpu: ['STELLAR_CPU', 'gkeyll-ci-stellar_cpu', 'continuous-integration/jenkins/stellar_cpu'],
        perlmutter_gpu: ['PERLMUTTER_GPU', 'gkeyll-ci-perlmutter_gpu', 'continuous-integration/jenkins/perlmutter_gpu'],
        team: ['TEAM_WORKSTATION', 'gkeyll-ci-team-workstation', 'continuous-integration/jenkins/team-workstation']
    ]
    final Map<Long, Map> states = new ConcurrentHashMap<>()
    final Map<Long, Long> progressBase = new ConcurrentHashMap<>()
    final File directory = new File(Jenkins.get().rootDir, 'gkeyll-queue-status')
    final AtomicBoolean checking = new AtomicBoolean(false)
    final AtomicBoolean refreshingPositions = new AtomicBoolean(false)
    final AtomicBoolean positionRefreshScheduled = new AtomicBoolean(false)
    final Set<Long> positionWrites = ConcurrentHashMap.newKeySet()
    final List<Object> statusLocks = (0..<32).collect { new Object() }

    Object statusLock(long id) { statusLocks[(id % statusLocks.size()) as int] }

    boolean updatingPosition(Queue.Item item) {
        def run = item.task instanceof Job ? null : nodeRun(item)
        positionWrites.contains(run ? run.queueId : item.id)
    }
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
        timer.scheduleWithFixedDelay({ sweep() } as Runnable, 60, 60, TimeUnit.SECONDS)
    }

    void sweep() {
        // Executor threads otherwise run anonymously and cannot discover jobs
        // on controllers with anonymous read access disabled.
        def context = ACL.as2(ACL.SYSTEM2)
        try {
            checkSuperseded()
            refreshPositions()
            refreshRunning()
        } catch (Exception failure) {
            // An exception escaping a scheduled task suppresses future sweeps.
            LOG.warning("Cannot refresh CI reporting (${failure.class.simpleName}); will retry.")
        } finally {
            context.close()
        }
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

    static String queuedDescription(long id, int position) {
        int ahead = Math.max(0, position - 1)
        def waiting = ahead ? "waiting behind ${ahead} other job${ahead == 1 ? '' : 's'}" :
            'no jobs ahead; waiting for an executor'
        "Gkeyll CI queued: ${waiting}" + runSuffix(id)
    }

    static boolean ownsQueuedStatus(Map status, def id) {
        status?.state == 'pending' && (status.description == "Gkeyll CI queued (Jenkins queue #${id})." ||
            status.description ==~ /Gkeyll CI queued: position [1-9][0-9]* of [1-9][0-9]* \(Jenkins queue #${id}\)\./ ||
            (status.description?.startsWith('Gkeyll CI queued: ') && ownsStatus(status, id)))
    }

    static boolean ownsRunningStatus(Map status, def id) {
        status?.state == 'pending' && status.description?.startsWith('Gkeyll CI running: ') &&
            ownsStatus(status, id)
    }

    static boolean ownsStatus(Map status, def id) {
        // Accept older reporters during rolling upgrades. The submission ID is
        // monotonic across jobs; a build number alone cannot identify its owner.
        status?.description?.endsWith(" (Jenkins queue #${id}).") ||
            status?.description ==~ /.* \(Jenkins (?:build #[0-9]+; )?ID ${id}\)\./
    }

    static String runSuffix(long id, def buildNumber = null) {
        " (Jenkins ${buildNumber ? 'build #' + buildNumber + '; ' : ''}ID ${id})."
    }

    Map currentStatus(Job job, Map config, Map metadata) {
        for (int page = 1; ; page++) {
            def path = "commits/${metadata.CI_QUEUE_COMMIT}/status?per_page=100" + (page > 1 ? "&page=${page}" : '')
            def statuses = request(job, config, 'GET', path).statuses
            def current = statuses.find { it.context == metadata.CI_QUEUE_CONTEXT }
            if (current || statuses.size() < 100) return current
        }
    }

    static String minutes(long milliseconds) {
        long value = Math.max(0L, milliseconds).intdiv(60000L)
        value < 60 ? "${value} m" : "${value.intdiv(60)} hr ${value % 60} m"
    }

    static String runningDescription(long id, long elapsed, long estimated, String stage = '', def buildNumber = null) {
        elapsed = Math.max(0L, elapsed)
        String progress = estimated <= 0 ? 'ETA unavailable' :
            (elapsed >= estimated ? 'ETA unavailable (estimate exceeded)' :
                "ETA ~${minutes(estimated - elapsed)}")
        String suffix = runSuffix(id, buildNumber)
        def timing = "elapsed ${minutes(elapsed)}; ${progress}"
        int room = Math.max(0, 140 - 'Gkeyll CI running: '.size() - timing.size() - suffix.size() - 2)
        def activity = stage && room ? stage.take(room) + '; ' : ''
        "Gkeyll CI running: ${activity}${timing}".take(140 - suffix.size()) + suffix
    }

    static String runningStage(Map status) {
        if (!status?.description?.startsWith('Gkeyll CI running: ')) return ''
        def stage = status.description.substring('Gkeyll CI running: '.size())
            .replaceFirst(/ \(Jenkins .*\)\.$/, '')
            .replaceFirst(/; (?:elapsed|started) .*/, '').replaceFirst(/\.$/, '')
        stage.startsWith('elapsed ') || stage.startsWith('started ') ? '' : stage
    }

    void refreshRunning() {
        try {
            Jenkins.get().getAllItems(Job).each { job ->
                if (!configuration(job)) return
                // Jenkins' in-progress chain avoids loading the full build history.
                for (def run = job.lastBuild; run != null; run = run.previousBuildInProgress) {
                    if (run.isBuilding()) refreshRun(run)
                }
            }
        } catch (Exception failure) {
            LOG.warning("Cannot refresh running CI (${failure.class.simpleName}); will retry.")
        }
    }

    void refreshRun(Run run) {
        synchronized (statusLock(run.queueId)) {
            try {
                def metadata = parameters(run)
                if (!run.isBuilding() || metadata.CI_QUEUE_STARTED != 'true' ||
                    metadata.CI_QUEUE_ID != run.queueId.toString() ||
                    !(metadata.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/)) return
                def config = configuration(run.parent)
                if (!config?.credential) return
                def current = currentStatus(run.parent, config, metadata)
                if (current && !ownsQueuedStatus(current, run.queueId) && !ownsRunningStatus(current, run.queueId)) return
                def environment = run.getAction(org.jenkinsci.plugins.workflow.cps.EnvActionImpl)?.getOverriddenEnvironment() ?: [:]
                // Command detail survives minute updates, including Pipelines
                // already running with an older reporting adapter.
                def stage = environment.CI_PROGRESS_STAGE ?: runningStage(current) ?: environment.CI_FAILURE_STAGE ?: 'starting'
                def description = runningDescription(run.queueId,
                    System.currentTimeMillis() - run.startTimeInMillis, run.estimatedDuration, stage, run.number)
                if (!run.isBuilding() || current?.description == description) return
                if (current?.id != null) progressBase.putIfAbsent(run.queueId, current.id as long)
                def body = [state: 'pending', context: metadata.CI_QUEUE_CONTEXT, description: description]
                if (current?.target_url) body.target_url = current.target_url
                request(run.parent, config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}", body)
            } catch (Exception failure) {
                LOG.warning("Queue #${run.queueId}: progress update failed (${failure.class.simpleName}); will retry.")
            }
        }
    }

    // Rank by original submission ID within each platform. A Pipeline's first
    // node() wait retains its original rank instead of becoming a new job.
    Map queuePositions() {
        def waiting = []
        Jenkins.get().queue.items.each { item ->
            if (item instanceof Queue.BuildableItem && item.isPending()) return // Already assigned an executor.
            def run = item.task instanceof Job ? null : nodeRun(item)
            def job = run?.parent ?: (item.task instanceof Job ? item.task : null)
            def config = job ? configuration(job) : null
            if (!config) return
            if (run) {
                def values = parameters(run)
                if (values.CI_QUEUE_ID != run.queueId.toString() || values.CI_QUEUE_STARTED == 'true') return
            }
            waiting << [id: run ? run.queueId : item.id, itemId: item.id, job: job, run: run, config: config]
        }
        def positions = [:]
        waiting.groupBy { it.config.platform }.values().each { entries ->
            def ordered = entries.unique { it.id }.sort { it.id }
            ordered.eachWithIndex { entry, index ->
                positions[entry.id] = entry + [position: index + 1, total: ordered.size()]
            }
        }
        positions
    }

    void schedulePositions() {
        if (timer.isShutdown() || !positionRefreshScheduled.compareAndSet(false, true)) return
        // Coalesce a burst of queue events, and keep HTTP requests off queue callbacks.
        timer.schedule({
            positionRefreshScheduled.set(false)
            refreshPositions()
        } as Runnable, 1, TimeUnit.SECONDS)
    }

    void refreshPositions() {
        if (!refreshingPositions.compareAndSet(false, true)) return
        try {
            queuePositions().values().each { entry ->
                synchronized (statusLock(entry.id)) {
                    try {
                        // Preparation publishes first, before releasing the job to Jenkins.
                        if (!entry.run && !states.get(entry.id)?.ready?.get()) return
                        def metadata = entry.run ? parameters(entry.run) : queuedParameters(entry.id, entry.job)
                        if (!(metadata.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/)) return
                        def current = currentStatus(entry.job, entry.config, metadata)
                        if (current && !ownsQueuedStatus(current, entry.id)) return
                        // Recompute after the HTTP call: the build may now be running or cancelled.
                        String description = null
                        Queue.withLock({
                            def live = queuePositions()[entry.id]
                            if (!live) return
                            def updated = queuedDescription(entry.id, live.position)
                            if (current?.description == updated) return
                            description = updated
                            positionWrites.add(entry.id)
                        } as Runnable)
                        if (description == null) return
                        try {
                            request(entry.job, entry.config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}",
                                [state: 'pending', context: metadata.CI_QUEUE_CONTEXT, description: description])
                        } finally {
                            positionWrites.remove(entry.id)
                            Jenkins.get().queue.scheduleMaintenance()
                        }
                    } catch (Exception failure) {
                        LOG.warning("Queue #${entry.id}: position update failed (${failure.class.simpleName}); will retry.")
                    }
                }
            }
        } catch (Exception failure) {
            LOG.warning("Cannot refresh queue positions (${failure.class.simpleName}); will retry.")
        } finally {
            refreshingPositions.set(false)
        }
    }

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
        def builder = HttpRequest.newBuilder(URI.create("https://api.github.com/repos/gkeyllorg/gkeyll/${path}"))
            .timeout(Duration.ofSeconds(25)).header('Accept', 'application/vnd.github+json')
            .header('Authorization', "Bearer ${credential.password.plainText}")
            .header('X-GitHub-Api-Version', '2022-11-28').header('User-Agent', 'gkeyll-jenkins-reporting')
            .header('Content-Type', 'application/json')
        def publisher = body == null ? HttpRequest.BodyPublishers.noBody() :
            HttpRequest.BodyPublishers.ofString(JsonOutput.toJson(body))
        def response = HTTP.send(builder.method(method, publisher).build(), HttpResponse.BodyHandlers.ofString())
        if (response.statusCode() < 200 || response.statusCode() >= 300) throw new IOException("GitHub HTTP ${response.statusCode()}")
        response.body() ? new JsonSlurper().parseText(response.body()) : [:]
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
            schedulePositions()
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
            CI_QUEUE_STARTED: 'false', CI_QUEUE_CANCEL_DESCRIPTION: '',
            CI_QUEUE_ENQUEUED_MS: saved.CI_QUEUE_ENQUEUED_MS ?: item.inQueueSince.toString(),
            CI_AGENT_STARTED_MS: '', CANDIDATE_PR: selected.pr,
            CANDIDATE_REF: selected.ref, BASELINE_REF: selected.baseline]
        // Do not change queue ParametersAction: doing so breaks Jenkins' duplicate
        // submission folding (e.g. every multibranch scan could enqueue another run).
        saveParameters(item.id, item.task, metadata)
        record.metadata = metadata
        if (record.cancelled.get()) return
        // onEnterWaiting fires before Jenkins publishes its new queue snapshot.
        // Wait for that short scheduling transaction, without holding the lock for HTTP.
        def position = null
        Queue.withLock({ position = queuePositions()[item.id] } as Runnable)
        if (record.cancelled.get()) return
        // Some queue callbacks run before the new snapshot is visible even
        // after withLock. Always publish the initial accepted state; the next
        // position sweep replaces this provisional rank once the item appears.
        position = position ?: [position: 1, total: 1]
        request(item.task, config, 'POST', "statuses/${selected.sha}",
            [state: 'pending', context: metadata.CI_QUEUE_CONTEXT,
             description: queuedDescription(item.id, position.position)])
    }

    void left(Queue.LeftItem item) {
        if (!item.cancelled) return
        def record = states.get(item.id)
        if (record != null) record.cancelled.set(true)
        workers.submit({
            // Wait for any in-flight pending POST before sending cancellation.
            synchronized (record ?: new Object()) {
                def metadata = record?.metadata ?: queuedParameters(item.id, item.task)
                finish(item.task, metadata, item.id, 'cancelled',
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
            def update = new ParametersAction([new StringParameterValue('CI_QUEUE_STARTED', 'true'),
                new StringParameterValue('CI_AGENT_STARTED_MS', System.currentTimeMillis().toString())], INTERNAL)
            run.addOrReplaceAction(run.getAction(ParametersAction).merge(update))
            workers.submit({ run.save(); refreshRun(run) } as Runnable)
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

    static String durationText(long milliseconds) {
        long seconds = Math.max(0L, milliseconds).intdiv(1000L)
        seconds < 60 ? "${seconds} s" : "${seconds.intdiv(60)} m ${seconds % 60} s"
    }

    static String escapeReport(String value) {
        (value ?: '').replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
    }

    static String reportExcerpt(String value, int byteLimit) {
        // Bound rendered UTF-8, including HTML expansion, not just input chars.
        def text = value ?: ''
        while (escapeReport(text).getBytes('UTF-8').length > byteLimit) {
            text = text.take(text.length().intdiv(2))
        }
        escapeReport(text) + (text == (value ?: '') ? '' : '\n[excerpt truncated]')
    }

    Map completionDetails(Run run, Map metadata) {
        def environment = run.getAction(org.jenkinsci.plugins.workflow.cps.EnvActionImpl)?.getOverriddenEnvironment() ?: [:]
        long end = run.startTimeInMillis + run.duration
        long start = (metadata.CI_AGENT_STARTED_MS ?: environment.CI_EXECUTION_START_MS ?: run.startTimeInMillis).toString().toLong()
        long queued = metadata.CI_QUEUE_ENQUEUED_MS ? Math.max(0L, start - metadata.CI_QUEUE_ENQUEUED_MS.toLong()) : -1L
        def stage = environment.CI_FAILURE_STAGE ?: 'Pipeline startup'
        def log = 'Jenkins console excerpt unavailable.'
        try { log = run.getLog(100).join('\n') }
        catch (Exception ignored) { /* Status publication must survive missing logs. */ }
        boolean bootstrap = environment.CI_BOOTSTRAP_COMPLETE != 'true' && environment.CI_BOOTSTRAP_COMMAND
        def excerpt = bootstrap ? (environment.CI_BOOTSTRAP_ERROR ?: log) : log
        def cause = environment.CI_FAILURE_MESSAGE ?: log.readLines().find { it.startsWith('ERROR: ') } ?: 'See diagnostic excerpt.'
        if (excerpt =~ /(?i)(unable to create thread|unable to create threaded lstat|Resource temporarily unavailable)/) {
            cause = 'Process/thread creation failed (resource exhaustion).'
        }
        [stage: stage, cause: cause, elapsed: durationText(end - start), queued: queued,
         started: start, ended: end, log: excerpt.take(12000), bootstrap: bootstrap,
         resources: bootstrap ? environment.CI_BOOTSTRAP_RESOURCES ?: 'Snapshot unavailable.' : '',
         command: bootstrap ? environment.CI_BOOTSTRAP_COMMAND : '',
         attempts: bootstrap ? environment.CI_BOOTSTRAP_ATTEMPT : '',
         exit: bootstrap ? environment.CI_BOOTSTRAP_EXIT : '',
         stageMs: bootstrap && environment.CI_BOOTSTRAP_STAGE_START_MS ?
            Math.max(0L, (environment.CI_BOOTSTRAP_STAGE_END_MS ?: end).toString().toLong() - environment.CI_BOOTSTRAP_STAGE_START_MS.toLong()) : null,
         retryMs: (environment.CI_BOOTSTRAP_RETRY_WAIT_MS ?: '0').toLong()]
    }

    String fallbackReport(Run run, Map config, Map metadata, Map detail) {
        // Per-run marker permits idempotent delivery without overwriting another
        // run's comment. No agent processes or workspace access are required.
        def marker = "<!-- gkeyll-ci-fallback ${metadata.CI_QUEUE_CONTEXT} ${run.queueId} -->"
        def body = "${marker}\n\n**Gkeyll CI ${run.result} after ${detail.elapsed}**\n\n" +
            "Jenkins build #${run.number}; ID ${run.queueId}. Candidate: `${metadata.CI_QUEUE_COMMIT}`.\n\n" +
            "**Stage:** ${reportExcerpt(detail.stage, 1000)}\n\n**Cause:** ${reportExcerpt(detail.cause, 2000)}\n\n" +
            "Execution: ${detail.elapsed}. Queue wait: ${detail.queued < 0 ? 'not recorded' : durationText(detail.queued)}. " +
            "Retry waiting: ${durationText(detail.retryMs)}.\n\n" +
            "Started: ${java.time.Instant.ofEpochMilli(detail.started)}; ended: ${java.time.Instant.ofEpochMilli(detail.ended)}.\n\n"
        if (detail.bootstrap) {
            body += "CI setup failed; compilation and tests did not start. Attempts: ${detail.attempts}; exit: ${detail.exit}. " +
                "Stage duration: ${detail.stageMs == null ? 'not recorded' : durationText(detail.stageMs)}.\n\n" +
                "**Command:**\n<pre>${reportExcerpt(detail.command, 6000)}</pre>\n\n" +
                "**Resource snapshot** (best effort; visible limits only):\n<pre>${reportExcerpt(detail.resources, 12000)}</pre>\n\n"
        }
        body += "**Diagnostic excerpt** (bounded; full log retained by Jenkins):\n<pre>${reportExcerpt(detail.log, 16000)}</pre>\n\n" +
            'The controller posted this fallback because no detailed Pipeline report was available.'
        new File(run.rootDir, 'gkeyll-fallback-report.md').setText(body, 'UTF-8')
        def path = "commits/${metadata.CI_QUEUE_COMMIT}/comments"
        def existing = null
        for (int page = 1; ; page++) {
            def comments = request(run.parent, config, 'GET', "${path}?per_page=100&page=${page}")
            existing = comments.find { it.body?.startsWith(marker + '\n') }
            if (existing || comments.size() < 100) break
        }
        def comment = existing ?: request(run.parent, config, 'POST', path, [body: body])
        def target = comment.html_url ?: ''
        new File(run.rootDir, 'gkeyll-fallback-delivery.json').setText(JsonOutput.toJson([comment: true, url: target]), 'UTF-8')
        target
    }

    void completed(Run run) {
        def metadata = parameters(run)
        if (metadata.CI_QUEUE_ID != run.queueId.toString()) return
        workers.submit({
            def state = run.result == Result.SUCCESS ? 'success' :
                (run.result in [Result.FAILURE, Result.UNSTABLE] ? 'failure' : 'cancelled')
            if (run.getAction(jenkins.model.InterruptedBuildAction)?.causes?.any { it.class.simpleName == 'ExceededTimeout' }) state = 'timed_out'
            def detail = completionDetails(run, metadata)
            def outcome = state == 'success' ? 'Passed' : state == 'timed_out' ? 'timed out' : state == 'cancelled' ? 'Cancelled' : 'Failed'
            def description = metadata.CI_QUEUE_CANCEL_DESCRIPTION ?:
                (state == 'cancelled' && metadata.CI_QUEUE_STARTED != 'true'
                    ? 'Gkeyll CI cancelled while queued.' : state == 'success' ? "Passed after ${detail.elapsed}." : "${outcome} after ${detail.elapsed}: ${detail.stage}; ${detail.cause}")
            finish(run.parent, metadata, run.queueId, state, description, run.number, run)
            progressBase.remove(run.queueId)
        } as Runnable)
    }

    void finish(Job job, Map metadata, long id, String state, String description, def buildNumber = null, Run run = null) {
        // Cancellation must follow an in-flight position POST, never precede it.
        synchronized (statusLock(id)) {
            finishStatus(job, metadata, id, state in ['cancelled', 'timed_out'] ? 'error' : state, description, buildNumber, run)
        }
    }

    void finishStatus(Job job, Map metadata, long id, String state, String description, def buildNumber = null, Run run = null) {
        if (metadata.CI_QUEUE_ID != id.toString() || !(metadata.CI_QUEUE_COMMIT ==~ /[0-9a-f]{40}/)) return
        def config = configuration(job)
        if (!config?.credential) return
        try {
            def current = currentStatus(job, config, metadata)
            def suffix = buildNumber ? runSuffix(id, buildNumber) : " (Jenkins queue #${id})."
            def failedAfterReport = "Jenkins failed after the test report: ${description}".take(140 - suffix.size()) + suffix
            // Artifact archival/cleanup can fail after tests were reported green.
            // Correct only this run's success, retaining the diagnostic report link.
            if (current?.state == 'success' && state != 'success' &&
                    ownsStatus(current, id)) {
                request(job, config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}",
                    [state: state, context: metadata.CI_QUEUE_CONTEXT, target_url: current.target_url,
                     description: failedAfterReport])
                return
            }
            // Preserve the Pipeline's detailed terminal status and another run's newer status.
            if (current && !ownsQueuedStatus(current, id) && !ownsRunningStatus(current, id)) return
            // The Pipeline can publish its detailed result while a progress POST
            // is in flight. Recover that result rather than leaving stale pending
            // or replacing its diagnostics with our generic completion message.
            if (ownsRunningStatus(current, id) && progressBase.containsKey(id)) {
                def history = request(job, config, 'GET', "commits/${metadata.CI_QUEUE_COMMIT}/statuses?per_page=100")
                def terminal = history.find { it.context == metadata.CI_QUEUE_CONTEXT &&
                    (it.id as long) > progressBase[id] && it.state in ['success', 'failure', 'error'] }
                if (terminal) {
                    def body = [state: terminal.state, context: metadata.CI_QUEUE_CONTEXT,
                        description: terminal.description]
                    if (terminal.state == 'success' && state != 'success' && ownsStatus(terminal, id)) {
                        body.state = state
                        body.description = failedAfterReport
                    }
                    if (terminal.target_url) body.target_url = terminal.target_url
                    request(job, config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}", body)
                    return
                }
            }
            def target = ''
            if (run && state != 'success') {
                try {
                    target = fallbackReport(run, config, metadata, completionDetails(run, metadata))
                } catch (Exception failure) {
                    new File(run.rootDir, 'gkeyll-fallback-delivery.json').setText(
                        JsonOutput.toJson([comment: false, error: failure.class.simpleName]), 'UTF-8')
                    LOG.warning("Queue #${id}: fallback comment failed (${failure.class.simpleName}); posting status anyway.")
                    description += '; report unavailable'
                }
            }
            // Comment delivery can take seconds. Recheck ownership after it so
            // an intervening detailed report or newer run keeps its status.
            if (run && state != 'success') {
                def latest = currentStatus(job, config, metadata)
                if (latest && !ownsQueuedStatus(latest, id) && !ownsRunningStatus(latest, id)) return
            }
            def body = [state: state, context: metadata.CI_QUEUE_CONTEXT,
                description: description.take(140 - suffix.size()) + suffix]
            if (target) body.target_url = target
            request(job, config, 'POST', "statuses/${metadata.CI_QUEUE_COMMIT}", body)
        } catch (Exception failure) {
            LOG.warning("Queue #${id}: GitHub final queue status failed (${failure.class.simpleName}).")
        }
    }
}

class GkeyllQueueListener extends QueueListener {
    final GkeyllQueueStatus service
    GkeyllQueueListener(GkeyllQueueStatus service) { this.service = service }
    @Override void onEnterWaiting(Queue.WaitingItem item) { service.begin(item); service.schedulePositions() }
    @Override void onLeft(Queue.LeftItem item) {
        if (item.task instanceof Job) {
            if (service.configuration(item.task)) service.left(item)
        } else {
            def run = GkeyllQueueStatus.nodeRun(item)
            if (run && service.configuration(run.parent)) service.nodeLeft(item, run)
        }
        service.schedulePositions()
    }
}

class GkeyllQueueDispatcher extends QueueTaskDispatcher {
    final GkeyllQueueStatus service
    GkeyllQueueDispatcher(GkeyllQueueStatus service) { this.service = service }
    @Override CauseOfBlockage canRun(Queue.Item item) {
        def record = service.begin(item)
        (record != null && !record.ready.get()) || service.updatingPosition(item) ? new CauseOfBlockage() {
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
