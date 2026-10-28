// Run with groovy.ui.GroovyMain, like test_stellar_resources.groovy.
// Exercise the trusted helper and GNU Make's effective architecture flags.
def source = new File('ci/jenkins/jenkinsfile.stellar-ai-cpu').text
def helpers = ['ciArchitectureFlags', 'applyCiArchitecture'].collect { name ->
    def start = source.indexOf("def ${name}(")
    assert start >= 0
    source.substring(start, source.indexOf('\n}\n', start) + 3)
}.join('\n')
def directory = File.createTempDir('stellar-ai-architecture-', '')
def environment = [:]
def run = { List command ->
    def builder = new ProcessBuilder(command.collect { it.toString() })
        .directory(directory).redirectErrorStream(true)
    builder.environment().putAll(environment)
    def process = builder.start()
    def output = process.inputStream.text
    assert process.waitFor() == 0 : output
    output.trim()
}
def pipeline = new GroovyShell(new Binding([
    withEnv: { List values, Closure body ->
        values.each { value ->
            def fields = value.toString().split('=', 2)
            environment[fields[0]] = fields[1]
        }
        body()
    },
    sh: { String script -> run(['bash', '-c', script]) }
])).parse(helpers)
try {
    def original = new File(directory, 'source-config.mak')
    def generated = new File(directory, 'config.mak')
    new File(directory, 'Makefile').text = '''include config.mak
.PHONY: flags
flags:
	@printf '%s\\n' '$(ARCH_FLAGS)'
'''
    ['ARCH_FLAGS = -march=native', 'override ARCH_FLAGS := -march=graniterapids'].each { setting ->
        original.text = "${setting}\nPREFIX = /fixture\n"
        def unchanged = original.text
        generated.text = original.text
        pipeline.applyCiArchitecture()
        assert original.text == unchanged
        assert generated.text.contains('PREFIX = /fixture')
        assert run(['make', '--no-print-directory', 'flags']) == pipeline.ciArchitectureFlags()
        assert run(['make', '--no-print-directory', 'flags', 'ARCH_FLAGS=-march=native']) == pipeline.ciArchitectureFlags()
    }
    println('Stellar AI CI architecture tests passed')
} finally {
    directory.deleteDir()
}
