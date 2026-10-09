# Jenkins CI tests

[CI overview](../README.md)

Run these commands from the repository root. The controller integration tests
require disposable Jenkins homes; never install their hooks on a production
controller.

## Reporting and queue listener

Run offline reporter tests with:

```sh
python3 -m unittest discover -s ci/jenkins/tests -p 'test_*.py'
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_jenkins_reporting.groovy
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_bootstrap_reporting.groovy
```

The queue integration test mocks GitHub for all four platforms, checking
supersession before agent allocation, cancellation, active runs, persistence,
reporting failures, queue counts, and progress without anonymous job access.
It needs Credentials, Credentials Binding, Folders, Pipeline: Job, Pipeline:
Groovy, Pipeline: Basic Steps, and Pipeline: Nodes and Processes, including
dependencies. With a local Jenkins WAR and plugin archives:

```sh
queue_test_home=$(mktemp -d)
mkdir -p "$queue_test_home/plugins" "$queue_test_home/init.groovy.d"
cp /path/to/test-plugin-archives/*.jpi "$queue_test_home/plugins/"
cp ci/jenkins/tests/test_queue_status.groovy "$queue_test_home/init.groovy.d/90-queue-test.groovy"
JENKINS_HOME="$queue_test_home" java \
  -Djenkins.install.runSetupWizard=false -Dgkeyll.queue.test=true \
  -Dgkeyll.queue.source="$PWD/ci/jenkins" \
  -jar /path/to/jenkins.war --httpListenAddress=127.0.0.1 --httpPort=18089
cat "$queue_test_home/queue-test-result.txt"
```

The test requires an empty job directory and shuts down its JVM on completion.

## Prebuilt dependencies

The offline config tests use Python 3, Groovy 2.4, GNU Make, and a C compiler:

```sh
java -cp /path/to/groovy-all.jar groovy.ui.GroovyMain ci/jenkins/tests/test_prebuilt_config.groovy
python3 -m unittest -v ci.jenkins.tests.test_prebuilt_config
```

The sandbox integration test uses a fresh Jenkins home with Pipeline: Job,
Pipeline: Groovy, Pipeline: Basic Steps, and Pipeline: Nodes and Processes
(including dependencies). With the WAR and plugin archives available locally:

```sh
prebuilt_test_home=$(mktemp -d)
mkdir -p "$prebuilt_test_home/plugins" "$prebuilt_test_home/init.groovy.d"
cp /path/to/test-plugin-archives/*.jpi "$prebuilt_test_home/plugins/"
cp ci/jenkins/tests/test_prebuilt_config_sandbox.groovy "$prebuilt_test_home/init.groovy.d/90-prebuilt-test.groovy"
JENKINS_HOME="$prebuilt_test_home" java \
  -Djenkins.install.runSetupWizard=false -Dgkeyll.prebuilt.test=true \
  -Dgkeyll.prebuilt.source="$PWD/ci/jenkins" \
  -jar /path/to/jenkins.war --httpPort=-1
cat "$prebuilt_test_home/prebuilt-test-result.txt"
```

This runs every helper in the default Groovy sandbox for candidate and baseline,
with and without MPI, then shuts down its JVM. Offline Groovy tests alone do
not check Jenkins sandbox permissions.

The installed-tool fixture checks configure/load, installation isolation, and
Lua failure exit codes without running simulations:

```sh
GKEYLL=/path/to/gkeyll/bin/gkeyll python3 -m unittest -v ci.jenkins.tests.test_regression_config
```

## Regression scheduling and output

The scheduler fixtures use tiny shell/make jobs, not plasma simulations:

```sh
LUAJIT=/path/to/luajit python3 -m unittest discover -s ci/jenkins/tests -p 'test_*.py'
```

The BGK output fixture also exercises the real CBC input for one step in serial
and on four MPI ranks, and compares the source rate and equilibrium arrays.
Set `GKEYLL` to the installed executable, `GKEYLL_CBC` to the built
`rt_gk_cbc_3x2v_p1`, and `MPIEXEC` to the matching MPI launcher, then run
`python3 -m unittest ci.jenkins.tests.test_bgk_source_io`. When using build-tree
executables, expose their shared libraries through `LD_LIBRARY_PATH`.

