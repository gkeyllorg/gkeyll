// Shared trusted Pipeline reporting. Load before the candidate checkout and
// keep the Python source in memory so the candidate cannot replace the publisher.
import groovy.transform.Field
import groovy.json.JsonOutput

@Field Map settings = [:]
@Field String pythonSource = ''
@Field String stageHistory = ''

def configure(Map values) {
    settings = values + [workspace: env.WORKSPACE]
    pythonSource = readFile('github_report.py')
    if (settings.reportingCommit) echo "CI reporting tools commit: ${settings.reportingCommit}"
    return this
}

def runReporter(String arguments, boolean authenticated = false) {
    // A Pipeline temporary directory is outside the candidate checkout. Python
    // isolated mode excludes candidate modules and PYTHONPATH from imports.
    def script = "${pwd(tmp: true)}/github_report.py"
    writeFile file: script, text: pythonSource
    withEnv(["CI_REPORT_SCRIPT=${script}", "CI_REPORT_CONTEXT=${settings.context}",
             "CI_REPORT_PLATFORM=${settings.platform}", "CI_REPORT_PR=${settings.pr ?: ''}",
             "CI_REPORT_REF=${settings.ref ?: ''}", "CI_REPORT_COMMIT=${settings.commit ?: ''}",
             "CI_REPORT_START_MS=${currentBuild.startTimeInMillis}",
             "CI_REPORT_END_MS=${new Date().time}",
             "CI_QUEUE_ID=${params.CI_QUEUE_ID ?: ''}"]) {
        def command = 'python3 -I "$CI_REPORT_SCRIPT" ' + arguments
        if (authenticated) {
            withCredentials([usernamePassword(credentialsId: settings.credential,
                    usernameVariable: 'GITHUB_USERNAME', passwordVariable: 'GITHUB_TOKEN')]) {
                sh command
            }
        } else {
            sh command
        }
    }
}

def publish(String result, String description = '', boolean detailed = false) {
    // Commands can execute in cached candidate/baseline trees; report artifacts
    // always belong to the build workspace, never to one of those shared trees.
    dir(settings.workspace) {
        // Replaced by the Python publisher when it runs. Retain an honest
        // delivery record even if credential binding or report generation fails.
        writeFile file: 'ci-report-delivery.json', text: JsonOutput.toJson([
            commit: settings.commit ?: '', context: settings.context, result: result,
            comment: false, status: false,
            errors: ['GitHub publication did not complete; see the Jenkins build log.']])
        withEnv(["CI_REPORT_RESULT=${result}", "CI_REPORT_DESCRIPTION=${description}"]) {
            if (detailed) {
                writeFile file: 'ci-stage-history.txt', text: stageHistory
                // Recreate provenance after candidate checkout deletes the workspace.
                if (env.CI_TRUSTED_CI_COMMIT?.trim()) {
                    writeFile file: 'ci-trusted-ci-commit.txt', text: "${env.CI_TRUSTED_CI_COMMIT.trim()}\n"
                }
                if (settings.reportingCommit) {
                    writeFile file: 'ci-reporting-commit.txt', text: "${settings.reportingCommit}\n"
                }
                runReporter('build --platform "$CI_REPORT_PLATFORM" --context "$CI_REPORT_CONTEXT" ' +
                    '--result "$CI_REPORT_RESULT" --pr "$CI_REPORT_PR" --output ci-report.md')
                description = readFile('ci-status-description.txt').trim()
            }
            // Keep a local report even if checkout failed before resolving a SHA.
            if (!(settings.commit ==~ /[0-9a-f]{40}/)) {
                echo 'WARNING: CI report retained locally; the candidate SHA is unavailable.'
                return
            }
            withEnv(["CI_REPORT_DESCRIPTION=${description}"]) {
                runReporter('update --platform "$CI_REPORT_PLATFORM" --context "$CI_REPORT_CONTEXT" ' +
                    '--result "$CI_REPORT_RESULT" --commit "$CI_REPORT_COMMIT" --pr "$CI_REPORT_PR" ' +
                    '--description "$CI_REPORT_DESCRIPTION"', true)
            }
        }
    }
    return description
}

def progress(String description) {
    // The controller polls this independently of the agent's long shell steps.
    env.CI_PROGRESS_STAGE = description.replaceFirst(/^Running: /, '').replaceFirst(/\.$/, '')
    try { publish('pending', description) }
    catch (err) {
        if (err.toString().contains('FlowInterruptedException')) throw err
        echo "WARNING: CI progress could not be published (${err.class.simpleName}); the build continues."
    }
}

def ciStage(String name, Closure body) {
    env.CI_FAILURE_STAGE = name
    // Keep history in Pipeline state as checkout deletes workspace files.
    stageHistory += "${new Date().format("yyyy-MM-dd'T'HH:mm:ss'Z'", TimeZone.getTimeZone('UTC'))}  ${name}\n"
    dir(settings.workspace) {
        writeFile file: 'ci-stage-history.txt', text: stageHistory
    }
    if (settings.commit) {
        progress("Running: ${name}.")
    }
    stage(name) { body.call() }
}

def writeCiFailureSummary(String result, String stageName, String message) {
    def compact = (message ?: 'unknown failure').replaceAll('[\\r\\n]+', ' ').trim().take(500)
    writeFile file: 'ci-failure-summary.txt', text: "result=${result}\nstage=${stageName ?: 'unknown'}\nmessage=${compact}\n"
}

def failureResult(def failure, String currentResult) {
    if (failure.toString().contains('FlowInterruptedException')) {
        if (failure.causes.any { it.class.simpleName == 'ExceededTimeout' }) return 'timed_out'
        return 'cancelled'
    }
    return currentResult == 'ABORTED' ? 'cancelled' : 'failure'
}

// Capture stdout/stderr and retain the original exit code. Never call these
// helpers with credential-bearing commands.
def loggedSh(String label, String script) {
    progress("Running: ${env.CI_FAILURE_STAGE ?: label} (${label}).")
    def quoted = "'" + script.replace("'", "'\"'\"'") + "'"
    withEnv(["CI_COMMAND_LOG=${env.WORKSPACE}/ci-command-logs/${label}.log"]) {
        sh """#!/bin/bash
            set -uo pipefail
            mkdir -p "\$WORKSPACE/ci-command-logs"
            bash ${script.startsWith('#!/bin/bash -l') ? '-l' : ''} -e -o pipefail -c ${quoted} 2>&1 | tee "\$CI_COMMAND_LOG"
            result=\$?
            printf '%s\\n' "\$result" > "\$CI_COMMAND_LOG.exit"
            exit "\$result"
        """
    }
}

def timeCommand(String metricFile, String command) {
    def logFile = metricFile.replaceFirst(/-seconds\.txt$/, '') + '.log'
    def detailFile = metricFile.replaceFirst(/[^\/]*$/, '') + 'ci-failure-detail.txt'
    def quoted = "'" + command.replace("'", "'\"'\"'") + "'"
    progress('Running: ' + metricFile.tokenize('/').last().replace('-seconds.txt', '').replace('-', ' ') + '.')
    sh """#!/bin/bash
        set -uo pipefail
        started=\$(date +%s)
        bash -e -o pipefail -c ${quoted} 2>&1 | tee '${logFile}'
        result=\$?
        printf '%s\\n' "\$result" > '${logFile}.exit'
        if [ "\$result" -ne 0 ]; then
            { printf 'Command: %s\\nFull log: %s (archived with the build)\\n' ${quoted} '${logFile}'; tail -n 100 '${logFile}'; } > '${detailFile}'
            exit "\$result"
        fi
        printf '%s\\n' "\$(( \$(date +%s) - started ))" > '${metricFile}'
    """
}

return this
