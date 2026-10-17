// Offline fault injection for bootstrap retries. No Jenkins or GitHub access.
class FlowInterruptedException extends RuntimeException { }
def source = new File('ci/jenkins/jenkinsfile.personal').text
source = source.substring(0, source.indexOf('// Manually selected local CI.')) + '\nreturn this\n'
def fixture = { List results, boolean interrupted = false ->
    def files = [:], environment = [:], sleeps = [], calls = []
    def binding = new Binding([
        env: environment,
        pwd: { Map options -> '/tmp/workspace@tmp' },
        writeFile: { Map value -> files[value.file.toString()] = value.text },
        readFile: { String path -> files[path] ?: '' },
        echo: { Object message -> },
        stage: { String name, Closure body -> body() },
        timeout: { Map options, Closure body -> assert options.time == 10; body() },
        sleep: { Map options -> sleeps << options.time; if (interrupted) throw new FlowInterruptedException() },
        error: { String message -> throw new IllegalStateException(message) },
        withEnv: { List values, Closure body ->
            def previous = new HashMap(environment)
            values.each { String value -> def pair = value.split('=', 2); environment[pair[0]] = pair[1] }
            try { body() } finally { values.each { String value -> def key = value.split('=', 2)[0]; environment[key] = previous[key] } }
        },
        sh: { Map options ->
            if (options.script.contains('Process limit:')) return 0
            calls << options.script
            def result = results.remove(0)
            files[environment.CI_BOOTSTRAP_COMMAND_LOG] = result.log
            result.code
        }
    ])
    [script: new GroovyShell(binding).evaluate(source), files: files, env: environment, sleeps: sleeps, calls: calls]
}
def resource = [code: 128, log: 'fatal: unable to create thread: Resource temporarily unavailable\nfatal: fetch-pack: invalid index-pack output\n']
def f = fixture([resource, resource, [code: 0, log: 'recovered']])
f.script.bootstrapStage('Load CI pipeline') { f.script.bootstrapSh('load-pipeline', 'git fetch') }
assert f.calls.size() == 3 && f.sleeps == [15, 30]
assert f.env.CI_FAILURE_MESSAGE == ''
assert f.files.findAll { k, v -> k.endsWith('.txt') }.size() == 3
assert f.files['/tmp/workspace@tmp/bootstrap-records/current/ci-bootstrap/load-pipeline-1.txt'].contains('exit=128')
assert f.env.CI_BOOTSTRAP_STAGE_END_MS.toLong() >= f.env.CI_BOOTSTRAP_STAGE_START_MS.toLong()
f = fixture([resource, resource, resource, resource])
try { f.script.bootstrapSh('load-pipeline', 'git fetch'); assert false }
catch (IllegalStateException expected) { assert expected.message.contains('process/thread') }
assert f.calls.size() == 4 && f.sleeps == [15, 30, 60]
f = fixture([[code: 128, log: "fatal: couldn't find remote ref missing"]])
try { f.script.bootstrapSh('load-pipeline', 'git fetch'); assert false }
catch (IllegalStateException expected) { assert expected.message.contains('128') }
assert f.calls.size() == 1 && f.sleeps.empty
f = fixture([resource], true)
try { f.script.bootstrapSh('load-pipeline', 'git fetch'); assert false }
catch (FlowInterruptedException expected) { }
assert f.calls.size() == 1
assert f.env.CI_BOOTSTRAP_RETRY_WAIT_MS.toLong() >= 0
println 'PASS: setup retries preserve attempts and timing, recover, exhaust, reject permanent failures and respect interruption'
