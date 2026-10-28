// Run with: java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_stellar_resources.groovy
// Exercise resource requests without submitting Slurm jobs.
['stellar_cpu', 'stellar-ai-cpu'].each { platform ->
def source = new File("ci/jenkins/jenkinsfile.${platform}").text
// Exercise the same scheduling invariants with each platform's environment.
source = source.replace('STELLAR_AI_CPU', 'STELLAR_CPU')
new GroovyShell().parse(source)
def helpers = ['buildJobs', 'positiveInteger', 'regressionJobs', 'regressionTestTimeout',
               'submitRegressionTests', 'submitParallelRegressionTests'].collect { name ->
    def start = source.indexOf("def ${name}(")
    assert start >= 0
    source.substring(start, source.indexOf('\n}\n', start) + 3)
}.join('\n')
def requests = []
def env = [WORKSPACE: '/fixture']
def binding = new Binding([
    env: env,
    error: { String message -> throw new IllegalArgumentException(message) },
    submitSlurmJob: { Object... args -> requests << args.toList() }
])
def pipeline = new GroovyShell(binding).parse(helpers)
assert pipeline.buildJobs() == '8'
assert pipeline.regressionJobs() == '16'
['4', '8', '16', '32'].each { budget ->
    env.STELLAR_CPU_REGRESSION_JOBS = budget
    ['baseline-create', 'candidate-check'].each { mode ->
        requests.clear()
        pipeline.submitRegressionTests(mode, '/baseline', '/candidate')
        pipeline.submitParallelRegressionTests(mode, '/baseline', '/candidate')
        def serial = requests[0], mpi = requests[1]
        assert serial[2] == budget : 'Serial allocation must match workers'
        assert serial[4].contains("CI_REGRESSION_JOBS=${budget}")
        assert mpi[2] == '1'
        assert mpi[5] == Math.max(4, budget.toInteger()).toString()
        assert mpi[4].contains("CI_REGRESSION_JOBS=${budget}")
    }
}
env.STELLAR_CPU_BUILD_JOBS = '6'
assert pipeline.buildJobs() == '6'
env.STELLAR_CPU_REGRESSION_JOBS = '0'
try {
    pipeline.regressionJobs()
    assert false : 'Zero worker budget must be rejected'
} catch (IllegalArgumentException expected) { }
println('Stellar resource request tests passed')

}
