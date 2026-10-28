// Offline compatibility tests for independently pinned reporting revisions.
['stellar_cpu', 'perlmutter_gpu', 'personal'].each { platform ->
    def source = new File("ci/jenkins/jenkinsfile.${platform}").text
    new GroovyShell().parse(source)
    def start = source.indexOf('def statisticsArchiveExcludes(')
    assert start >= 0
    def helper = source.substring(start, source.indexOf('\n}\n', start) + 3)
    def warnings = []
    def binding = new Binding([echo: { String message -> warnings << message }])
    def pipeline = new GroovyShell(binding).parse(helper)
    def archive = source.readLines().find { it.contains('excludes: statisticsArchiveExcludes()') }
    assert archive && archive.contains('archiveArtifacts')

    // Current helper returns exclusions only after successful consolidation.
    binding.setVariable('reporting', [statisticsArchiveExcludes: { '*-seconds.txt' }])
    assert pipeline.statisticsArchiveExcludes() == '*-seconds.txt'
    binding.setVariable('reporting', [statisticsArchiveExcludes: { '' }])
    assert pipeline.statisticsArchiveExcludes() == ''
    assert warnings.empty

    // Older loaded scripts fail differently in plain Groovy and Jenkins CPS.
    binding.setVariable('reporting', new GroovyShell().parse('return this'))
    assert pipeline.statisticsArchiveExcludes() == ''
    binding.setVariable('reporting', [statisticsArchiveExcludes: {
        throw new NoSuchMethodError("No such DSL method 'statisticsArchiveExcludes' found among steps")
    }])
    assert pipeline.statisticsArchiveExcludes() == ''
    assert warnings.size() == 2

    // Do not mask unrelated linkage errors, failed archival, or interruption.
    [new NoSuchMethodError('unrelated helper'),
     new MissingMethodException('unrelated', String, [] as Object[]),
     new RuntimeException('interrupted')].each { failure ->
        binding.setVariable('reporting', [statisticsArchiveExcludes: { throw failure }])
        try {
            pipeline.statisticsArchiveExcludes()
            assert false: 'Unrelated failure must propagate'
        } catch (Throwable actual) {
            assert actual.is(failure)
        }
    }
}
println('Statistics archive compatibility tests passed')
