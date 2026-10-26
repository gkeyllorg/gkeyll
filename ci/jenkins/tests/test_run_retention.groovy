// Run with groovy.ui.GroovyMain. No Jenkins controller or Slurm jobs required.
def directory = File.createTempDir('gkeyll-retention-', '')
try {
    def bin = new File(directory, 'bin')
    bin.mkdirs()
    def python = new File(bin, 'python3')
    python.text = '''#!/bin/sh
printf '%s\n' "$@" > "$CI_CAPTURE_ARGS"
'''
    python.setExecutable(true)
    def tools = new File(directory, 'trusted tools')
    tools.mkdirs()
    new File(tools, 'module_load.sh').text = ':\n'
    def captured = new File(directory, 'arguments')
    ['personal', 'stellar_cpu', 'perlmutter_gpu', 'stellar-ai-cpu'].each { platform ->
        def source = new File("ci/jenkins/jenkinsfile.${platform}").text
        new GroovyShell().parse(source)
        def name = platform == 'stellar-ai-cpu' ? 'pruneRunDirectories' : 'cacheManager'
        def start = source.indexOf("def ${name}(")
        assert start >= 0
        def helper = source.substring(start, source.indexOf('\n}\n', start) + 3)
        [null, '', '/custom/auth file'].each { auth ->
            def environment = [
                PATH: "${bin}:${System.getenv('PATH')}",
                CI_CAPTURE_ARGS: captured.path,
                GKEYLL_CI_ROOT: new File(directory, 'CI root').path,
                GKEYLL_CI_BUILDS_TO_KEEP: '3',
                CI_CACHE_PLATFORM: platform,
                CI_CACHE_MANAGER: new File(tools, 'cache_manager.py').path,
                CI_BASELINE_CACHE_TOOL: new File(tools, 'baseline_cache.sh').path,
                CI_STELLAR_AI_CPU_TOOLS_DIR: tools.path,
                BUILD_TAG: 'jenkins-test-18', JOB_NAME: 'test', BUILD_NUMBER: '18',
                BUILD_URL: 'http://browser.invalid/job/test/18/'
            ].collectEntries { key, value -> [(key): value.toString()] }
            if (auth != null) environment.CI_JENKINS_READ_AUTH_FILE = auth
            def pipeline = new GroovyShell(new Binding([
                env: environment,
                withEnv: { List values, Closure body ->
                    def saved = new LinkedHashMap(environment)
                    values.each { value ->
                        def fields = value.toString().split('=', 2)
                        environment[fields[0]] = fields[1]
                    }
                    try { body() } finally { environment.clear(); environment.putAll(saved) }
                },
                sh: { Object step ->
                    def script = step instanceof Map ? step.script : step
                    def builder = new ProcessBuilder('bash', '--noprofile', '--norc', '-c', script.toString())
                        .redirectErrorStream(true)
                    builder.environment().putAll(environment)
                    def process = builder.start()
                    def output = process.inputStream.text
                    assert process.waitFor() == 0 : output
                    output
                }
            ])).parse(helper)
            if (platform == 'stellar-ai-cpu') pipeline.pruneRunDirectories(environment.GKEYLL_CI_ROOT, '3')
            else pipeline.cacheManager('finish')
            def arguments = captured.text.split('\n', -1).toList()
            def option = { String flag ->
                assert arguments.contains(flag) : "Missing ${flag}: ${arguments}"
                arguments[arguments.indexOf(flag) + 1]
            }
            assert option('--runs-to-keep') == '3'
            assert option('--root') == environment.GKEYLL_CI_ROOT
            assert option('--tag') == environment.BUILD_TAG
            assert option('--retain-runs') == '0'
            def defaultAuth = platform == 'stellar-ai-cpu' ? "${environment.GKEYLL_CI_ROOT}/jenkins_home/jenkins-cli.auth" : ''
            assert option('--auth-file') == (auth == null ? defaultAuth : auth)
            assert arguments.contains('--run-layout')
        }
    }

    // Execute the actual Stellar AI finally block with failing external steps.
    def source = new File('ci/jenkins/jenkinsfile.stellar-ai-cpu').text
    def start = source.indexOf('            def slurmFinalized = false')
    def end = source.indexOf('\n        }\n        }\n    }\n}', start)
    assert start >= 0 && end > start
    def finalizer = source.substring(start, end)
    ['none', 'slurm', 'archive', 'prune'].each { failure ->
        def pruned = 0, cleaned = 0
        def binding = new Binding([
            env: [:], currentBuild: [:], terminalStatus: 'success', terminalDescription: '', failureMessage: '',
            candidateDir: '', baselineCommit: '', candidateCommit: '', candidatePr: '', candidateRef: '',
            gkeyllCiRoot: '/fixture', buildsToKeep: '2',
            finalizeUnfinishedSlurmJob: { if (failure == 'slurm') throw new RuntimeException('slurm') },
            writeCiTimingSummary: {}, timingStatusDescription: { status -> status },
            withEnv: { values, body -> body() }, regressionJobs: { '16' },
            writeCiFailureSummary: { Object... ignored -> }, echo: { message -> },
            statisticsArchiveExcludes: { '' },
            archiveArtifacts: { options -> if (failure == 'archive') throw new RuntimeException('archive') },
            pruneRunDirectories: { root, count ->
                assert root == '/fixture' && count == '2'
                pruned++
                if (failure == 'prune') throw new RuntimeException('prune')
            },
            deleteDir: { cleaned++ }
        ])
        try {
            new GroovyShell(binding).evaluate(finalizer)
            assert failure == 'none'
        } catch (RuntimeException error) {
            assert error.message == failure
        }
        assert pruned == (failure in ['none', 'prune'] ? 1 : 0)
        assert cleaned == 1
    }
    println('Pipeline run-retention tests passed')
} finally {
    directory.deleteDir()
}
