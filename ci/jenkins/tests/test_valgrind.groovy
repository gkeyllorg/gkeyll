// Run with Groovy 2.4; exercise the trusted Pipeline helpers without Jenkins.
def source = new File('ci/jenkins/jenkins_reporting.groovy').text
def helpers = source.substring(source.indexOf('def valgrindEnabled()'), source.indexOf('def runReporter('))
def calls = []
def binding = new Binding([
    env: [WORKSPACE: '/workspace', CI_FAILURE_STAGE: 'Build candidate'],
    settings: [platform: 'personal', valgrindScript: '/trusted/valgrind.py'],
    error: { String message -> throw new IllegalArgumentException(message) },
    withEnv: { List values, Closure body -> calls << values; body() },
    loggedSh: { String name, String command -> calls << [name, command] },
    timeCommand: { String metric, String command -> calls << [metric, command] }
])
binding.setVariable('ciStage', { String name, Closure body ->
    binding.env.CI_FAILURE_STAGE = name
    calls << name
    body()
})
def script = new GroovyShell(binding).parse(helpers)
assert !script.valgrindEnabled()
script.prepareValgrind('candidate')
script.runValgrind()
assert calls.empty
['2', 'true', '-1'].each { value ->
    binding.env.GKEYLL_USE_VALGRIND = value
    try { script.valgrindEnabled(); assert false }
    catch (IllegalArgumentException expected) { }
}
binding.env.GKEYLL_USE_VALGRIND = '1'
['personal', 'team', 'stellar_cpu'].each { platform ->
    binding.settings.platform = platform
    assert script.valgrindEnabled()
}
script.prepareValgrind('candidate')
script.prepareValgrind('baseline')
script.runValgrind()
assert calls.toString().contains('candidate-valgrind-configure')
assert calls.toString().contains('baseline-valgrind-configure')
assert calls.toString().contains('candidate-valgrind-seconds.txt')
assert calls.toString().contains('ci-valgrind/candidate')
assert binding.env.CI_FAILURE_STAGE == 'Build candidate'
binding.setVariable('timeCommand', { String metric, String command -> throw new IllegalStateException('fixture leak') })
try { script.runValgrind(); assert false }
catch (IllegalStateException expected) { }
assert binding.env.CI_FAILURE_STAGE == 'Valgrind serial unit tests'
binding.settings.platform = 'perlmutter_gpu'
try { script.valgrindEnabled(); assert false }
catch (IllegalArgumentException expected) { assert expected.message.contains('CPU') }
println 'PASS: Valgrind flag, CPU gating, configuration, failure stage, and artifact commands'
