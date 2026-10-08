// Run with: java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/test_prebuilt_config.groovy
// No Jenkins controller or dependency builds required.
import java.nio.file.Files

def root = new File('ci/jenkins')
def shell = new GroovyShell()
['personal', 'team_workstation', 'stellar_cpu', 'perlmutter_gpu'].each {
    shell.parse(new File(root, "jenkinsfile.${it}"))
}
def helpers = ['personal', 'stellar_cpu', 'perlmutter_gpu'].collect {
    def source = new File(root, "jenkinsfile.${it}").text
    source.substring(source.indexOf('def usePrebuiltConfig('), source.indexOf('\n}\n', source.indexOf('def usePrebuiltConfig(')) + 3)
}
assert helpers.toSet().size() == 1 : 'Keep trusted helper copies consistent'
def tmp = Files.createTempDirectory('gkeyll-prebuilt-test').toFile()
try {
    def original = new File(tmp, 'installed-config.mak')
    original.text = '''PREFIX=/opt/existing
INSTALL_PREFIX=/opt/existing
LAPACK_INC_DIR=$(PREFIX)/OpenBLAS/include
SUPERLU_LIB_DIR=${PREFIX}/superlu/lib
CONF_MPI_LIB_DIR=/opt/mpi/lib
CONF_LUA_LIB_DIR=/opt/PREFIX/lua/lib
'''
    def saved = original.text
    def output = new File(tmp, 'config.mak')
    def binding = new Binding([
        env: [:],
        fileExists: { String path -> new File(path).isFile() },
        readFile: { Map args -> new File(args.file).text },
        writeFile: { Map args -> new File(tmp, args.file).text = args.text },
        echo: { String message -> },
        error: { String message -> throw new IllegalArgumentException(message) }
    ])
    def helper = new GroovyShell(binding).parse(helpers[0])
    assert !helper.usePrebuiltConfig('/unused')
    assert !output.exists()
    ['relative/config.mak', '/nonexistent/gkeyll-config.mak'].each { path ->
        binding.env.GKEYLL_CI_PREBUILT_CONFIG = path
        try {
            helper.usePrebuiltConfig('/unused')
            assert false : 'Invalid config should fail'
        } catch (IllegalArgumentException expected) { }
    }
    binding.env.GKEYLL_CI_PREBUILT_CONFIG = original.absolutePath
    ['candidate', 'baseline'].each { tree ->
        def prefix = "${tmp}/${tree}/gkylsoft"
        assert helper.usePrebuiltConfig(prefix)
        def probe = new File(tmp, 'probe.mk')
        probe.text = '''include config.mak
.PHONY: probe
probe:
\t@printf '%s\\n' '$(PREFIX)' '$(INSTALL_PREFIX)' '$(LAPACK_INC_DIR)' '$(SUPERLU_LIB_DIR)' '$(CONF_MPI_LIB_DIR)' '$(CONF_LUA_LIB_DIR)'
'''
        def process = ['make', '-s', '-f', probe.name, 'probe'].execute(null, tmp)
        def result = process.text.readLines()
        assert process.waitFor() == 0
        assert result == [prefix, prefix, '/opt/existing/OpenBLAS/include', '/opt/existing/superlu/lib', '/opt/mpi/lib', '/opt/PREFIX/lua/lib']
        assert original.text == saved : 'Shared config must remain untouched'
    }
    // Exercise the actual local/team build function with Jenkins steps mocked.
    def personal = new File(root, 'jenkinsfile.personal').text
    def start = personal.indexOf('def buildTree(')
    def buildTree = personal.substring(start, personal.indexOf('\n}\n', start) + 3)
    def setupCalls = []
    def buildCalls = []
    binding.setVariable('withEnv', { List values, Closure body -> body.call() })
    binding.setVariable('loggedSh', { String label, String command -> setupCalls << command })
    binding.setVariable('timeCommand', { String metric, String command -> buildCalls << command })
    binding.env.WORKSPACE = tmp.absolutePath
    binding.env.PERSONAL_MPIEXEC = '/opt/mpi/bin/mpiexec'
    def pipeline = new GroovyShell(binding).parse(helpers[0] + buildTree)
    ['candidate', 'baseline'].each { tree ->
        pipeline.buildTree(tree, "${tmp}/${tree}/gkylsoft", '', '', '2')
    }
    assert setupCalls.empty : 'Prebuilt mode must skip both machine scripts'
    assert buildCalls.size() == 8 : 'Both trees must still build, install and test'
    binding.env.GKEYLL_CI_PREBUILT_CONFIG = ''
    pipeline.buildTree('candidate', "${tmp}/normal/gkylsoft", 'mkdeps.test.sh', 'configure.test.sh', '2')
    assert setupCalls.size() == 1
    assert setupCalls[0].contains('./machines/mkdeps.test.sh')
    assert setupCalls[0].contains('./machines/configure.test.sh')
    assert buildCalls.size() == 12
    binding.env.GKEYLL_CI_PREBUILT_CONFIG = original.absolutePath
    original.text = 'CC=cc\n'
    try {
        helper.usePrebuiltConfig('/unused')
        assert false : 'Config without PREFIX should fail'
    } catch (IllegalArgumentException expected) { }
    println 'PASS: Jenkinsfiles parse; prebuilt config validation, dependency paths and candidate/baseline isolation'
} finally {
    tmp.deleteDir()
}
