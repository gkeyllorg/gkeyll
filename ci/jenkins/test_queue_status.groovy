// Integration test: run ONLY in a disposable Jenkins home (see README.md).
// GitHub calls are mocked; no real credentials, statuses, or builds are used.
import com.cloudbees.hudson.plugins.folder.Folder
import hudson.ExtensionList
import hudson.init.InitMilestone
import hudson.model.ParametersAction
import hudson.model.Job
import hudson.model.ParametersDefinitionProperty
import hudson.model.Queue
import hudson.model.Result
import hudson.model.StringParameterDefinition
import hudson.model.StringParameterValue
import hudson.slaves.EnvironmentVariablesNodeProperty
import jenkins.model.Jenkins
import org.jenkinsci.plugins.workflow.cps.CpsFlowDefinition
import org.jenkinsci.plugins.workflow.job.WorkflowJob
import org.jenkinsci.plugins.workflow.job.properties.DisableConcurrentBuildsJobProperty
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.CountDownLatch
import java.util.concurrent.TimeUnit

if (System.getProperty('gkeyll.queue.test') != 'true') {
    throw new IllegalStateException('This test requires -Dgkeyll.queue.test=true and a disposable Jenkins home.')
}
def source = new File(System.getProperty('gkeyll.queue.source'))
def jenkins = Jenkins.get()
assert jenkins.allItems.empty : 'Use a fresh, disposable Jenkins home.'
def shell = new GroovyShell(jenkins.pluginManager.uberClassLoader)
def service = shell.evaluate(new File(source, 'queue-status.groovy'))
service.timer.shutdownNow() // Tests trigger deterministic sweeps.

def heads = new ConcurrentHashMap<String, String>()
def statuses = new ConcurrentHashMap<String, Map>()
def calls = Collections.synchronizedList([])
def history = Collections.synchronizedList([])
def statusIds = new java.util.concurrent.atomic.AtomicLong()
def controls = new ConcurrentHashMap()
service.metaClass.request = { Job job, Map config, String method, String path, Map body = null ->
    calls << [method: method, path: path, body: body]
    if (controls.fail) throw new IOException('Simulated GitHub outage')
    if (path == 'pulls/10' && controls.gate) {
        controls.entered.countDown()
        assert controls.gate.await(30, TimeUnit.SECONDS)
    }
    if (method == 'POST' && controls.postGate && path == 'statuses/' + 'e' * 40) {
        controls.postEntered.countDown()
        assert controls.postGate.await(30, TimeUnit.SECONDS)
    }
    if (method == 'POST' && controls.positionGate && path == 'statuses/' + '40' * 20) {
        controls.positionEntered.countDown()
        assert controls.positionGate.await(30, TimeUnit.SECONDS)
    }
    if (method == 'POST') {
        assert path.startsWith('statuses/')
        def sha = path.substring('statuses/'.length())
        if (controls.progressTerminal && body.description.startsWith('Gkeyll CI running:')) {
            history << [sha: sha, id: statusIds.incrementAndGet(), context: body.context,
                state: 'failure', description: 'Detailed failure during progress POST', target_url: 'https://example.test/build/1']
            controls.remove('progressTerminal')
        }
        def status = new HashMap(body) + [id: statusIds.incrementAndGet()]
        history << status + [sha: sha]
        statuses[sha + ':' + body.context] = status
        return [:]
    }
    if (path.endsWith('/statuses?per_page=100')) {
        return history.findAll { it.sha == path.split('/')[1] }.sort { -it.id }.take(100)
    }
    if (path.startsWith('pulls/')) return [head: [sha: heads[path.substring(6)]], base: [ref: 'main']]
    if (path.endsWith('/status?per_page=100')) {
        def sha = path.split('/')[1]
        return [statuses: statuses.findAll { k, v -> k.startsWith(sha + ':') }.values().toList()]
    }
    if (path.startsWith('commits/')) return [sha: controls.branchSha ?: 'd' * 40]
    throw new AssertionError("Unexpected request: ${method} ${path}")
}

Thread.start('gkeyll-queue-integration-tests') {
    int checks = 0
    def await = { String label, Closure condition ->
        long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(40)
        while (!condition()) {
            assert System.nanoTime() < deadline : "Timed out: ${label}"
            Thread.sleep(100)
        }
    }
    def passed = { String label -> checks++; println("PASS: ${label}") }
    try {
        assert service.runningDescription(42L, 29L * 60000, 139L * 60000) ==
            'Gkeyll CI running: started 29 m ago; ~20%; est. remaining 1 hr 50 m (Jenkins queue #42).'
        assert service.runningDescription(42L, 60000L, -1L).contains('ETA unavailable')
        assert service.runningDescription(42L, 60000L, 0L).contains('ETA unavailable')
        assert service.runningDescription(42L, 60000L, 60000L).contains('exceeded estimate')
        assert !service.runningDescription(42L, 120000L, 60000L).contains('100%')
        assert service.runningDescription(42L, -1L, 60000L).contains('started 0 m ago; ~0%')
        assert service.runningDescription(Long.MAX_VALUE, Long.MAX_VALUE, 1L).size() <= 140
        assert service.ownsRunningStatus([state: 'pending', description:
            service.runningDescription(Long.MAX_VALUE, Long.MAX_VALUE, 1L)], Long.MAX_VALUE)
        passed('progress formatting covers estimates, missing history, overruns, clock skew and the status length limit')
        await('Jenkins startup') { jenkins.initLevel == InitMilestone.COMPLETED }
        jenkins.setNumExecutors(1)
        jenkins.setLabelString('queue-test-runner')
        def vars = [PERSONAL_GITHUB_CREDENTIAL_ID: 'mock', PERSONAL_STATUS_CONTEXT: 'ci/personal-test',
            STELLAR_CPU_GITHUB_CREDENTIAL_ID: 'mock', PERLMUTTER_GPU_GITHUB_CREDENTIAL_ID: 'mock',
            TEAM_WORKSTATION_GITHUB_CREDENTIAL_ID: 'mock']
        jenkins.globalNodeProperties.add(new EnvironmentVariablesNodeProperty(
            vars.collect { k, v -> new EnvironmentVariablesNodeProperty.Entry(k, v) }))
        def create = { parent, String name, String script ->
            def job = parent.createProject(WorkflowJob, name)
            job.definition = new CpsFlowDefinition(script, false)
            job.addProperty(new ParametersDefinitionProperty(
                ['CANDIDATE_PR', 'CANDIDATE_REF', 'BASELINE_REF'].collect { new StringParameterDefinition(it, '') }))
            job.save()
            job
        }
        def submit = { job, int quiet, Map values ->
            def future = values ? job.scheduleBuild2(quiet, new ParametersAction(
                values.collect { k, v -> new StringParameterValue(k, v) })) : job.scheduleBuild2(quiet)
            assert future != null
            await('queue entry') { jenkins.queue.items.any { it.task == job } }
            jenkins.queue.items.findAll { it.task == job }.max { it.id }
        }
        def pending = { item, String sha, String context ->
            await("pending queue ${item.id}") {
                service.states[item.id]?.ready?.get() &&
                    service.ownsQueuedStatus(statuses[sha + ':' + context], item.id)
            }
        }
        def absent = { long id -> !(jenkins.queue.getItem(id) instanceof Queue.Item) || jenkins.queue.getItem(id) instanceof Queue.LeftItem }
        def contexts = ['ci/personal-test', 'continuous-integration/jenkins/stellar_cpu',
            'continuous-integration/jenkins/perlmutter_gpu', 'continuous-integration/jenkins/team-workstation']
        def team = jenkins.createProject(Folder, 'gkeyll-ci-team-workstation')
        def jobs = [create(jenkins, 'gkeyll-ci-personal', "echo 'BUILD_EXECUTED'"),
            create(jenkins, 'gkeyll-ci-stellar_cpu', "echo 'BUILD_EXECUTED'"),
            create(jenkins, 'gkeyll-ci-perlmutter_gpu', "echo 'BUILD_EXECUTED'"),
            create(team, 'PR-4', "echo 'BUILD_EXECUTED'")]
        def items = []
        jobs.eachWithIndex { job, i ->
            heads[(i + 1).toString()] = (i + 1).toString() * 40
            def item = submit(job, 600, i == 3 ? [:] : [CANDIDATE_PR: (i + 1).toString()])
            items << item
            pending(item, heads[(i + 1).toString()], contexts[i])
            assert job.lastBuild == null
            assert service.queuedParameters(item.id, job).CI_QUEUE_COMMIT == heads[(i + 1).toString()]
        }
        passed('all four platforms publish pending before the Pipeline starts, including automatic PR discovery')
        assert submit(jobs[3], 600, [:]).id == items[3].id
        assert jenkins.queue.items.size() == 4
        passed('repeated automatic submissions still coalesce into one Jenkins queue item')

        assert statuses['1' * 40 + ':' + contexts[0]].description.contains('position 1 of 1')
        heads['40'] = '40' * 20
        def laterJob = create(team, 'PR-40', "echo 'BUILD_EXECUTED'")
        def later = submit(laterJob, 600, [:])
        pending(later, '40' * 20, contexts[3])
        service.refreshPositions()
        assert statuses['4' * 40 + ':' + contexts[3]].description.contains('position 1 of 2')
        assert statuses['40' * 20 + ':' + contexts[3]].description.contains('position 2 of 2')
        assert statuses['1' * 40 + ':' + contexts[0]].description.contains('position 1 of 1')
        def posts = calls.count { it.method == 'POST' }
        service.refreshPositions()
        assert calls.count { it.method == 'POST' } == posts
        passed('positions count waiting builds per platform and unchanged positions are not reposted')

        heads['1'] = 'a' * 40
        service.checkSuperseded()
        await('superseded queue status') { statuses['1' * 40 + ':' + contexts[0]]?.state == 'error' }
        assert absent(items[0].id)
        assert jobs[0].lastBuild == null
        assert statuses['1' * 40 + ':' + contexts[0]].description.contains('cancelled while queued; superseded')
        assert !absent(items[1].id)
        passed('a new PR head cancels the older queued build without running it or affecting another PR')

        assert jenkins.queue.cancel(jenkins.queue.getItem(items[1].id))
        await('manual cancellation') { statuses['2' * 40 + ':' + contexts[1]]?.state == 'error' }
        assert jobs[1].lastBuild == null
        passed('manual queue cancellation clears pending')

        def duplicate = submit(jobs[2], 600, [CANDIDATE_PR: '3', RETRY_REQUEST: 'second'])
        pending(duplicate, '3' * 40, contexts[2])
        service.checkSuperseded()
        assert !absent(items[2].id) && !absent(duplicate.id)
        assert jenkins.queue.cancel(jenkins.queue.getItem(items[2].id))
        await('old cancellation callback') { !service.states.containsKey(items[2].id) }
        assert service.ownsQueuedStatus(statuses['3' * 40 + ':' + contexts[2]], duplicate.id)
        assert jenkins.queue.cancel(jenkins.queue.getItem(duplicate.id))
        await('duplicate cancellation') { statuses['3' * 40 + ':' + contexts[2]]?.state == 'error' }
        passed('same-commit reruns are retained; cancelling an older queue item preserves the newer status')

        heads['4'] = 'b' * 40
        service.checkSuperseded()
        await('automatic PR cancellation') { statuses['4' * 40 + ':' + contexts[3]]?.state == 'error' }
        assert jobs[3].lastBuild == null
        passed('automatic PR builds are superseded without a replacement submission')
        controls.positionGate = new CountDownLatch(1)
        controls.positionEntered = new CountDownLatch(1)
        def refreshing = Thread.start { service.refreshPositions() }
        assert controls.positionEntered.await(20, TimeUnit.SECONDS)
        assert service.positionWrites.contains(later.id)
        def dispatcher = ExtensionList.lookup(hudson.model.queue.QueueTaskDispatcher)
            .find { it.class.simpleName == 'GkeyllQueueDispatcher' }
        assert dispatcher.canRun(jenkins.queue.getItem(later.id)) != null
        assert calls.last().body.description.contains('position 1 of 1')
        assert jenkins.queue.cancel(jenkins.queue.getItem(later.id))
        controls.positionGate.countDown()
        refreshing.join(30000)
        assert !refreshing.alive
        controls.remove('positionGate')
        await('positioned job cancellation') { statuses['40' * 20 + ':' + contexts[3]]?.state == 'error' }
        assert !service.positionWrites.contains(later.id)
        passed('positions advance after cancellation; an in-flight update cannot overwrite a later cancellation')

        def branch = submit(jobs[0], 600, [CANDIDATE_REF: 'feature/test', BASELINE_REF: 'main'])
        pending(branch, 'd' * 40, contexts[0])
        assert calls.any { it.path == 'commits/feature%2Ftest' }
        jenkins.queue.save()
        service.states.clear()
        controls.branchSha = 'e' * 40
        jenkins.queue.load()
        pending(branch, 'd' * 40, contexts[0])
        assert service.queuedParameters(branch.id, jobs[0]).CI_QUEUE_COMMIT == 'd' * 40
        assert jenkins.queue.cancel(jenkins.queue.getItem(branch.id))
        passed('branch selectors are URL-encoded and remain pinned after the queue is saved and reloaded')

        // Hold a GitHub request and prove that neither the build nor Jenkins' queue lock is held by it.
        heads['10'] = 'e' * 40
        controls.gate = new CountDownLatch(1)
        controls.entered = new CountDownLatch(1)
        def delayed = submit(jobs[0], 0, [CANDIDATE_PR: '10'])
        assert controls.entered.await(20, TimeUnit.SECONDS)
        assert jobs[0].lastBuild == null
        def unrelated = create(jenkins, 'unrelated', "echo 'unrelated build'")
        assert unrelated.scheduleBuild2(0).get(30, TimeUnit.SECONDS).result == Result.SUCCESS
        assert jenkins.queue.cancel(jenkins.queue.getItem(delayed.id))
        controls.gate.countDown()
        controls.remove('gate')
        await('cancelled preparation') { !service.states.containsKey(delayed.id) }
        passed('GitHub latency does not block Jenkins; cancellation during candidate resolution prevents execution')

        controls.postGate = new CountDownLatch(1)
        controls.postEntered = new CountDownLatch(1)
        def posting = submit(jobs[0], 600, [CANDIDATE_PR: '10'])
        assert controls.postEntered.await(20, TimeUnit.SECONDS)
        assert jenkins.queue.cancel(jenkins.queue.getItem(posting.id))
        controls.postGate.countDown()
        await('cancellation after in-flight pending POST') { statuses['e' * 40 + ':' + contexts[0]]?.state == 'error' }
        controls.remove('postGate')
        passed('cancellation cannot be overwritten by its delayed pending notification')

        heads['5'] = '5' * 40
        def waiting = create(team, 'PR-5', "node('offline-queue-test-agent') { echo 'BUILD_EXECUTED' }")
        waiting.scheduleBuild2(0)
        await('Pipeline waiting for its first node') {
            waiting.lastBuild?.isBuilding() && jenkins.queue.items.any { service.nodeRun(it) == waiting.lastBuild }
        }
        def waitRun = waiting.lastBuild
        assert statuses['5' * 40 + ':' + contexts[3]].state == 'pending'
        service.refreshPositions()
        assert service.queuePositions()[waitRun.queueId].position == 1
        assert statuses['5' * 40 + ':' + contexts[3]].description.contains('position 1 of 1')
        service.refreshRunning()
        assert service.ownsQueuedStatus(statuses['5' * 40 + ':' + contexts[3]], waitRun.queueId)
        heads['5'] = 'f' * 40
        service.checkSuperseded()
        await('waiting Pipeline aborted') { !waitRun.isBuilding() }
        await('waiting Pipeline cancellation status') { statuses['5' * 40 + ':' + contexts[3]]?.state == 'error' }
        assert waitRun.result == Result.ABORTED
        assert !waitRun.getLog(200).any { it == 'BUILD_EXECUTED' }
        assert statuses['5' * 40 + ':' + contexts[3]].description.contains('superseded')
        passed('supersession also cancels a Pipeline waiting for its first agent, before its node body executes')

        heads['6'] = '6' * 40
        def running = create(team, 'PR-6', "node('queue-test-runner') { echo 'BUILD_EXECUTED'; sleep time: 3600, unit: 'SECONDS' }")
        running.scheduleBuild2(0)
        await('node allocated') { running.lastBuild?.getLog(200)?.any { it == 'BUILD_EXECUTED' } }
        assert service.parameters(running.lastBuild).CI_QUEUE_STARTED == 'true'
        await('initial running status') { service.ownsRunningStatus(statuses['6' * 40 + ':' + contexts[3]], running.lastBuild.queueId) }
        service.refreshRunning()
        assert statuses['6' * 40 + ':' + contexts[3]].description.contains('ETA unavailable')
        posts = calls.count { it.method == 'POST' }
        service.refreshRunning()
        assert calls.count { it.method == 'POST' } == posts
        def runningKey = '6' * 40 + ':' + contexts[3]
        def savedProgress = statuses[runningKey]
        statuses[runningKey] = [state: 'success', context: contexts[3], description: 'Detailed terminal status']
        service.refreshRunning()
        assert statuses[runningKey].description == 'Detailed terminal status'
        statuses[runningKey] = [state: 'pending', context: contexts[3],
            description: service.queuedDescription(running.lastBuild.queueId + 1, 1, 1)]
        service.refreshRunning()
        assert service.ownsQueuedStatus(statuses[runningKey], running.lastBuild.queueId + 1)
        statuses[runningKey] = savedProgress + [description: "Gkeyll CI running: starting (Jenkins queue #${running.lastBuild.queueId})."]
        controls.progressTerminal = true
        service.refreshRunning()
        assert !controls.progressTerminal
        passed('running sweeps skip unchanged, completed and newer-run statuses; agent waits remain queued')
        heads['6'] = 'c' * 40
        service.checkSuperseded()
        assert running.lastBuild.isBuilding()
        running.lastBuild.doStop()
        await('running test cleanup') { !running.lastBuild.isBuilding() }
        await('detailed result restored after progress race') {
            statuses[runningKey].description == 'Detailed failure during progress POST'
        }
        assert statuses[runningKey].target_url == 'https://example.test/build/1'
        passed('a newer PR commit does not abort CI that has already acquired its agent')

        heads['11'] = 'a1' * 20
        heads['12'] = 'a2' * 20
        jobs[0].addProperty(new DisableConcurrentBuildsJobProperty())
        jobs[0].definition = new CpsFlowDefinition(
            "node('queue-test-runner') { echo params.CI_QUEUE_COMMIT; sleep time: 3600, unit: 'SECONDS' }", false)
        jobs[0].scheduleBuild2(0, new ParametersAction(new StringParameterValue('CANDIDATE_PR', '11')))
        await('serialized job active') { jobs[0].lastBuild?.getLog(200)?.any { it == 'a1' * 20 } }
        def active = jobs[0].lastBuild
        def blocked = submit(jobs[0], 0, [CANDIDATE_PR: '12'])
        pending(blocked, 'a2' * 20, contexts[0])
        assert jobs[0].lastBuild == active
        heads['12'] = 'a3' * 20
        service.checkSuperseded()
        await('serialized superseded status') { statuses['a2' * 20 + ':' + contexts[0]]?.state == 'error' }
        assert absent(blocked.id)
        assert active.isBuilding() && jobs[0].lastBuild == active
        active.doStop()
        await('serialized test cleanup') { !active.isBuilding() }
        passed('disableConcurrentBuilds does not delay queue reporting or supersession; pinned SHA reaches Pipeline params')

        heads['7'] = '7' * 40
        def failed = create(team, 'PR-7', "error 'intentional pre-agent failure'")
        assert failed.scheduleBuild2(0).get(30, TimeUnit.SECONDS).result == Result.FAILURE
        await('pre-agent failure status') { statuses['7' * 40 + ':' + contexts[3]]?.state == 'failure' }
        passed('a failure before agent allocation closes the pending status')

        heads['8'] = '8' * 40
        def completed = create(team, 'PR-8', "sleep time: 2, unit: 'SECONDS'")
        def completion = completed.scheduleBuild2(0)
        await('pending before detailed result') { statuses['8' * 40 + ':' + contexts[3]]?.state == 'pending' }
        statuses['8' * 40 + ':' + contexts[3]] = [state: 'success', context: contexts[3], description: 'Detailed timings preserved']
        assert completion.get(30, TimeUnit.SECONDS).result == Result.SUCCESS
        assert statuses['8' * 40 + ':' + contexts[3]].description == 'Detailed timings preserved'
        passed('normal Pipeline result descriptions are preserved')

        controls.fail = true
        heads['9'] = '9' * 40
        def outage = create(team, 'PR-9', "echo 'reporter outage does not block builds'")
        assert outage.scheduleBuild2(0).get(30, TimeUnit.SECONDS).result == Result.SUCCESS
        controls.remove('fail')
        passed('GitHub outages release the dispatcher and allow the existing Pipeline to run')

        ['personal', 'team_workstation', 'stellar_cpu', 'perlmutter_gpu'].each {
            shell.classLoader.parseClass(new File(source, "jenkinsfile.${it}"))
        }
        passed('all four Jenkinsfiles compile with Jenkins Groovy')
        new File(jenkins.rootDir, 'queue-test-result.txt').text = "PASS: ${checks} integration checks\n"
        println("PASS: ${checks} integration checks")
        jenkins.cleanUp()
        System.exit(0)
    } catch (Throwable failure) {
        failure.printStackTrace()
        new File(jenkins.rootDir, 'queue-test-result.txt').text = "FAIL: ${failure}\n"
        System.exit(1)
    }
}
