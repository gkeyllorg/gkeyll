// Integration test: run ONLY in a disposable Jenkins home (see README.md).
// Exercise the trusted helper with real Pipeline steps and the default sandbox.
import hudson.init.InitMilestone
import hudson.model.Result
import hudson.slaves.EnvironmentVariablesNodeProperty
import jenkins.model.Jenkins
import org.jenkinsci.plugins.scriptsecurity.scripts.ScriptApproval
import org.jenkinsci.plugins.workflow.cps.CpsFlowDefinition
import org.jenkinsci.plugins.workflow.job.WorkflowJob
import java.util.concurrent.TimeUnit

if (System.getProperty('gkeyll.prebuilt.test') != 'true') {
    throw new IllegalStateException('This test requires -Dgkeyll.prebuilt.test=true and a disposable Jenkins home.')
}
def source = new File(System.getProperty('gkeyll.prebuilt.source'))
def jenkins = Jenkins.get()
assert jenkins.allItems.empty : 'Use a fresh, disposable Jenkins home.'
assert ScriptApproval.get().approvedSignatures.length == 0 : 'Use the default sandbox whitelist.'

Thread.start('gkeyll-prebuilt-sandbox-tests') {
    int status = 1
    def result = new File(jenkins.rootDir, 'prebuilt-test-result.txt')
    try {
        long deadline = System.nanoTime() + TimeUnit.SECONDS.toNanos(60)
        while (jenkins.initLevel != InitMilestone.COMPLETED) {
            assert System.nanoTime() < deadline : 'Timed out waiting for Jenkins startup.'
            Thread.sleep(100)
        }
        jenkins.setNumExecutors(1)
        jenkins.globalNodeProperties.add(new EnvironmentVariablesNodeProperty(
            new EnvironmentVariablesNodeProperty.Entry('CI_PREBUILT_TEST_SCRIPT',
                new File(source, 'prebuilt_config.py').absolutePath),
            new EnvironmentVariablesNodeProperty.Entry('CI_CACHE_TEST_SCRIPT',
                new File(source, 'baseline_cache.sh').absolutePath)))
        ['personal', 'stellar_cpu', 'perlmutter_gpu'].each { platform ->
            def pipeline = new File(source, "jenkinsfile.${platform}").text
            def start = pipeline.indexOf('def usePrebuiltConfig(')
            def helper = pipeline.substring(start, pipeline.indexOf('\n}\n', start) + 3)
            def reporting = new File(source, 'jenkins_reporting.groovy').text
            def envStart = reporting.indexOf('def withBaselineEnvironment(')
            helper += reporting.substring(envStart, reporting.indexOf('\n}\n', envStart) + 3)
            def cacheStart = reporting.indexOf('def withBaselineCache(')
            helper += reporting.substring(cacheStart, reporting.indexOf('\n}\n', cacheStart) + 3)
            [true, false].each { mpi ->
                def job = jenkins.createProject(WorkflowJob, "prebuilt-${platform}-${mpi ? 'mpi' : 'serial'}")
                job.definition = new CpsFlowDefinition('import groovy.transform.Field\n@Field def reporting\n@Field def settings\n'
                    + helper + "\ndef mpi = ${mpi}\n" + '''
node {
    reporting = [settings: [prebuiltScript: env.CI_PREBUILT_TEST_SCRIPT]]
    env.GKEYLL_CI_PREBUILT_CONFIG = ''
    env.MPI_HOME = ''
    env.OPAL_PREFIX = ''
    env.PERSONAL_MPI_HOME = ''
    env.PERSONAL_MPIEXEC = ''
    assert !usePrebuiltConfig('/unused')
    def root = pwd()
    settings = [cacheScript: env.CI_CACHE_TEST_SCRIPT]
    withBaselineCache("${root}/cache", 'personal') {
        assert fileExists("${root}/cache/baseline-cache/personal/.lock/owner")
    }
    assert !fileExists("${root}/cache/baseline-cache/personal/.lock")
    try {
        withBaselineCache("${root}/cache", 'personal') { error('lease cleanup fixture') }
        assert false
    } catch (Exception expected) { }
    assert !fileExists("${root}/cache/baseline-cache/personal/.lock")
    env.CI_RUN_DIR = "${root}/run"
    env.GKEYLL_CI_PREBUILT_CONFIG = "${root}/original.mak"
    env.LD_LIBRARY_PATH = '/existing/lib=with space'
    withEnv(["CI_TEST_MPI=${mpi ? '1' : '0'}"]) {
        sh \'\'\'#!/bin/bash
            set -euo pipefail
            mkdir -p original/OpenBLAS/{include,lib} original/openmpi/{include,lib,bin}
            mkdir -p original/gkeyll/share/adas
            printf 'ADAS fixture' > original/gkeyll/share/adas/ioniz_h.npy
            printf '#!/bin/sh\\nexit 0\\n' > original/openmpi/bin/mpiexec
            cat > original.mak <<EOF
PREFIX=$PWD/original
BUILD_APP=pkpm
LAPACK_INC_DIR=$PWD/original/OpenBLAS/include
LAPACK_LIB_DIR=$PWD/original/OpenBLAS/lib
USE_MPI=$CI_TEST_MPI
CONF_MPI_INC_DIR=$PWD/original/openmpi/include
CONF_MPI_LIB_DIR=$PWD/original/openmpi/lib
EOF
        \'\'\'
    }
    ['candidate', 'baseline'].each { tree ->
        def prepareTree = {
        dir(tree) {
            def prefix = "${root}/${tree}/gkylsoft"
            assert usePrebuiltConfig(prefix)
            assert readFile("${prefix}/gkeyll/share/adas/ioniz_h.npy") == 'ADAS fixture'
            assert readFile('config.mak').contains("PREFIX=${prefix}\\n")
            assert fileExists('alltargets.mak')
            def dependencies = tree == 'baseline' ? "${root}/baseline-cache/dependencies" : "${root}/run/dependencies"
            assert env.GKEYLL_CI_DEPENDENCY_ENV == "${dependencies}/env.sh"
            assert fileExists(env.GKEYLL_CI_DEPENDENCY_ENV)
            assert env.LD_LIBRARY_PATH.startsWith("${dependencies}/")
            assert env.LD_LIBRARY_PATH.endsWith(':/existing/lib=with space')
            if (mpi) {
                assert env.MPI_HOME == "${dependencies}/gkylsoft/openmpi"
                assert env.PERSONAL_MPI_HOME == env.MPI_HOME
                assert env.OPAL_PREFIX == env.MPI_HOME
                assert env.PERSONAL_MPIEXEC == "${env.MPI_HOME}/bin/mpiexec"
                assert fileExists(env.PERSONAL_MPIEXEC)
            } else {
                assert !env.MPI_HOME && !env.PERSONAL_MPIEXEC
            }
        }
        }
        if (tree == 'baseline') {
            def candidateEnvironment = env.GKEYLL_CI_DEPENDENCY_ENV
            withBaselineEnvironment("${root}/baseline-cache") { prepareTree() }
            assert env.GKEYLL_CI_DEPENDENCY_ENV == candidateEnvironment
            assert env.LD_LIBRARY_PATH.startsWith("${root}/run/dependencies/")
            try {
                withBaselineEnvironment("${root}/baseline-cache") {
                    env.LD_LIBRARY_PATH = '/baseline-only'
                    error('fixture interruption')
                }
                assert false
            } catch (Exception expected) { }
            assert env.LD_LIBRARY_PATH.startsWith("${root}/run/dependencies/")
        } else { prepareTree() }
    }
}
''', true)
                def run = job.scheduleBuild2(0).get(90, TimeUnit.SECONDS)
                assert run.result == Result.SUCCESS : run.getLog(200).join('\n')
                assert ScriptApproval.get().pendingSignatures.empty : 'No script approvals should be needed.'
                println "PASS: ${job.name} candidate and baseline in the default Groovy sandbox"
            }
        }
        result.text = 'PASS: all prebuilt helpers run in the default Groovy sandbox with and without MPI\n'
        status = 0
    } catch (Throwable failure) {
        result.text = "FAIL: ${failure}\n"
        failure.printStackTrace()
    } finally {
        System.exit(status)
    }
}
