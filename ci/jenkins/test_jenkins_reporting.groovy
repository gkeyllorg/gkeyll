// Offline tests for the trusted Pipeline adapter; no credentials or GitHub calls.
def files = ['github_report.py': '# trusted source', 'ci-status-description.txt': 'Passed.\n']
def commands = []
def credentials = []
def environments = []
def directories = []
def environment = [WORKSPACE: '/workspace', BUILD_NUMBER: '10', CI_TRUSTED_CI_COMMIT: 'b' * 40]
def binding = new Binding([
    env: environment, params: [CI_QUEUE_ID: '42'],
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
assert script.publish('success', '', true) == 'Passed.'
assert commands.any { it.contains(' build --platform ') }
assert commands.last().contains(' update --platform ')
assert commands.every { it.startsWith('python3 -I "$CI_REPORT_SCRIPT" ') }
assert directories.every { it == '/workspace' }
assert credentials.every { it.credentialsId == 'existing-token' }
assert environments.any { it == 'CI_QUEUE_ID=42' }
assert files['/workspace@tmp/github_report.py'] == '# trusted source'
assert files['ci-trusted-ci-commit.txt'] == 'b' * 40 + '\n'
assert files['ci-reporting-commit.txt'] == 'c' * 40 + '\n'

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
