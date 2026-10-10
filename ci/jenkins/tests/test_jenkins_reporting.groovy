// Offline tests for the trusted Pipeline adapter; no credentials or GitHub calls.
def files = ['github_report.py': '# trusted source', 'ci-status-description.txt': 'Passed.\n']
def commands = []
def credentials = []
def environments = []
def directories = []
def environment = [WORKSPACE: '/workspace', BUILD_NUMBER: '10', CI_TRUSTED_CI_COMMIT: 'b' * 40]
def binding = new Binding([
    env: environment, params: [CI_QUEUE_ID: '42'],
    currentBuild: [startTimeInMillis: 1700000000000L],
    readFile: { String path -> files[path] },
    writeFile: { Map value -> files[value.file] = value.text },
    fileExists: { String path -> files.containsKey(path) },
    pwd: { Map options -> '/workspace@tmp' },
    dir: { String path, Closure body -> directories << path; body() },
    withEnv: { List values, Closure body -> environments.addAll(values); body() },
    usernamePassword: { Map value -> value },
    withCredentials: { List values, Closure body -> credentials.addAll(values); body() },
    sh: { String command -> commands << command },
    stage: { String name, Closure body -> body() },
    echo: { String message -> }
])
def script = new GroovyShell(binding).evaluate(new File('ci/jenkins/jenkins_reporting.groovy'))
script.configure([platform: 'personal', context: 'ci/test', credential: 'existing-token',
    commit: 'a' * 40, reportingCommit: 'c' * 40])
script.ciStage('Build candidate') { }
assert files['ci-stage-history.txt'].contains('Build candidate')
assert environment.CI_FAILURE_STAGE == 'Build candidate'
assert environment.CI_PROGRESS_STAGE == 'Build candidate'
assert script.publish('success', '', true) == 'Passed.'
assert commands.any { it.contains(' build --platform ') }
assert commands.last().contains(' update --platform ')
assert commands.every { it.startsWith('python3 -I "$CI_REPORT_SCRIPT" ') }
assert directories.every { it == '/workspace' }
assert credentials.every { it.credentialsId == 'existing-token' }
assert environments.any { it == 'CI_QUEUE_ID=42' }
assert environments.any { it == 'CI_REPORT_START_MS=1700000000000' }
assert environments.any { it.startsWith('CI_REPORT_END_MS=') && it.split('=')[1].toLong() >= 1700000000000L }
assert files['/workspace@tmp/github_report.py'] == '# trusted source'
assert files['ci-trusted-ci-commit.txt'] == 'b' * 40 + '\n'
assert files['ci-reporting-commit.txt'] == 'c' * 40 + '\n'

try {
    script.ciStage('Failed stage') { throw new RuntimeException('stage fixture') }
    assert false
} catch (RuntimeException expected) { assert expected.message == 'stage fixture' }
assert script.stageTimings.last().result == 'failure'
assert script.stageTimings.last().elapsed_ms >= 0
assert environment.CI_FAILURE_MESSAGE.contains('stage fixture')
script.ciStage('Build candidate') { }

script.timeCommand('candidate-unit-build-seconds.txt', 'make unit')
assert environment.CI_PROGRESS_STAGE == 'candidate unit build'
assert environment.CI_FAILURE_STAGE == 'Build candidate'
script.loggedSh('configure-candidate', './configure')
assert environment.CI_PROGRESS_STAGE == 'Build candidate (configure-candidate)'

// Execute the generated wrapper: startup files may return nonzero, but the
// workload must still enforce errexit and pipefail and retain output/status.
def shellRoot = java.nio.file.Files.createTempDirectory('gkeyll-reporting-shell-').toFile()
def savedWorkspace = environment.WORKSPACE
try {
    environment.WORKSPACE = shellRoot.absolutePath
    def startup = new File(shellRoot, 'startup.sh')
    startup.text = 'false\n'
    [false, true].each { login ->
        [
            [name: 'success', body: "printf \"it's working\\n\"; echo stderr-message >&2", status: 0],
            [name: 'errexit', body: 'echo before-failure; false; echo should-not-run', status: 1],
            [name: 'pipefail', body: 'echo before-failure; false | cat; echo should-not-run', status: 1]
        ].each { fixture ->
            def label = "shell-${login}-${fixture.name}"
            def workload = (login ? '#!/bin/bash -l\n' : '#!/bin/bash\n') + fixture.body
            script.loggedSh(label, workload)
            def wrapper = new File(shellRoot, 'wrapper.sh')
            wrapper.text = commands.last()
            def processBuilder = new ProcessBuilder('/bin/bash', wrapper.absolutePath).redirectErrorStream(true)
            processBuilder.environment().putAll([
                WORKSPACE: shellRoot.absolutePath,
                CI_COMMAND_LOG: new File(shellRoot, "ci-command-logs/${label}.log").absolutePath,
                BASH_ENV: startup.absolutePath
            ])
            def process = processBuilder.start()
            def finished = process.waitFor(30, java.util.concurrent.TimeUnit.SECONDS)
            if (!finished) process.destroyForcibly()
            assert finished: 'shell fixture timed out'
            def output = process.inputStream.text
            assert process.exitValue() == fixture.status: output
            def log = new File(shellRoot, "ci-command-logs/${label}.log").text
            assert new File(shellRoot, "ci-command-logs/${label}.log.exit").text.trim() == fixture.status.toString()
            assert !log.contains('should-not-run')
            if (fixture.status == 0) {
                assert log.contains("it's working") && log.contains('stderr-message')
            } else {
                assert log.contains('before-failure')
            }
        }
    }
} finally {
    environment.WORKSPACE = savedWorkspace
    shellRoot.deleteDir()
}
println 'PASS: shell startup precedes strict mode; workload failures and logs are preserved'

// Generate and retain the failure report even when checkout never resolved a SHA.
script.settings.commit = ''
commands.clear()
credentials.clear()
script.publish('failure', 'Checkout failed', true)
assert commands.size() == 1 && commands.first().contains(' build --platform ')
assert credentials.empty
script.writeCiFailureSummary('failure', 'Build candidate', 'failed\non two lines')
assert files['ci-failure-summary.txt'].contains('message=failed on two lines')
assert script.failureResult(new RuntimeException('command failed'), 'FAILURE') == 'failure'
assert script.failureResult(new RuntimeException('aborted'), 'ABORTED') == 'cancelled'
println 'PASS: shared Pipeline publishing uses existing credentials and workspace artifacts, isolates Python, and retains early failure reports'

// Contending jobs wait; the lease is released even when a Pipeline body fails.
script.settings.cacheScript = '/trusted/baseline_cache.sh'
commands.clear()
def attempts = [75, 0]
binding.setVariable('timeout', { Map options, Closure body -> body() })
binding.setVariable('waitUntil', { Map options, Closure body -> while (!body()) { } })
binding.setVariable('error', { String message -> throw new IllegalStateException(message) })
binding.setVariable('sh', { Object command ->
    commands << command
    command instanceof Map ? attempts.remove(0) : null
})
try {
    script.withBaselineCache('/cache', 'personal') { throw new RuntimeException('cache body failed') }
    assert false
} catch (RuntimeException expected) { assert expected.message == 'cache body failed' }
assert attempts.empty
assert commands.size() == 3
assert commands.last().contains(' release ')
println 'PASS: cache contention waits and failed builds release the lease'
