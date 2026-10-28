// One-time controller bootstrap. Install as init.groovy.d/10-stellar-ai-cpu.groovy.
// Run only in a new, private Jenkins home; never resets an existing installation.
import jenkins.model.Jenkins
import jenkins.install.InstallState
import hudson.security.HudsonPrivateSecurityRealm
import hudson.security.FullControlOnceLoggedInAuthorizationStrategy
import hudson.slaves.EnvironmentVariablesNodeProperty
import hudson.model.*
import hudson.tasks.LogRotator
import org.jenkinsci.plugins.workflow.job.WorkflowJob
import org.jenkinsci.plugins.workflow.job.properties.DisableConcurrentBuildsJobProperty
import org.jenkinsci.plugins.workflow.cps.CpsScmFlowDefinition
import hudson.plugins.git.GitSCM
import hudson.plugins.git.BranchSpec
import hudson.plugins.git.UserRemoteConfig
import jenkins.security.ApiTokenProperty
import java.security.SecureRandom

final j = Jenkins.get()
final root = System.getenv('GKEYLL_CI_ROOT')
assert root?.startsWith('/') : 'GKEYLL_CI_ROOT must be an absolute shared path'
final marker = new File(j.rootDir, '.stellar-ai-cpu-configured')
if (marker.exists()) return
assert j.getAllItems().isEmpty() : 'Bootstrap requires a new controller'
assert !j.installState.isSetupComplete() : 'Bootstrap must not replace a configured controller'

// Write secrets only into the private controller home, never the console.
def secretFile = { String name, String value ->
    def path = new File(j.rootDir, name).toPath()
    java.nio.file.Files.writeString(path, value + '\n')
    java.nio.file.Files.setPosixFilePermissions(path, java.nio.file.attribute.PosixFilePermissions.fromString('rw-------'))
}
def passwordBytes = new byte[32]
new SecureRandom().nextBytes(passwordBytes)
def password = Base64.encoder.encodeToString(passwordBytes)
def realm = new HudsonPrivateSecurityRealm(false)
def user = realm.createAccount('admin', password)
j.securityRealm = realm
def strategy = new FullControlOnceLoggedInAuthorizationStrategy()
strategy.allowAnonymousRead = false
j.authorizationStrategy = strategy
j.setCrumbIssuer(new hudson.security.csrf.DefaultCrumbIssuer(false))
j.slaveAgentPort = -1
secretFile('admin-password', password)
def token = user.getProperty(ApiTokenProperty).tokenStore.generateNewToken('cluster-cli')
user.save()
secretFile('jenkins-cli.auth', 'admin:' + token.plainValue)

j.numExecutors = 1
j.mode = Node.Mode.EXCLUSIVE
j.labelString = 'stellar-ai-cpu'
def variables = new EnvironmentVariablesNodeProperty()
variables.envVars.putAll([
    GKEYLL_CI_ROOT: root,
    GKEYLL_CI_TRUSTED_REF: 'jenkins-github-more-reporting',
    GKEYLL_CI_BUILDS_TO_KEEP: '2',
    STELLAR_AI_CPU_GITHUB_CREDENTIAL_ID: 'gkeyll-github-stellar-ai-cpu',
    STELLAR_AI_CPU_SLURM_ACCOUNT: 'csd',
    STELLAR_AI_CPU_SLURM_QOS: 'short',
    STELLAR_AI_CPU_UNIT_QOS: 'test',
    STELLAR_AI_CPU_NODE_LABEL: 'stellar-ai-cpu',
    STELLAR_AI_CPU_BUILD_JOBS: '8',
    STELLAR_AI_CPU_REGRESSION_JOBS: '16'
])
j.globalNodeProperties.add(variables)
def job = j.createProject(WorkflowJob, 'gkeyll-ci-stellar-ai-cpu')
job.description = 'Manual Stellar AI CPU Slurm CI. Workflows always come from jenkins-github-more-reporting, independently of candidate/baseline refs.'
def scm = new GitSCM(
    [new UserRemoteConfig('https://github.com/gkeyllorg/gkeyll.git', null, null, null)],
    [new BranchSpec('*/jenkins-github-more-reporting')], false, [], null, null, [])
job.definition = new CpsScmFlowDefinition(scm, 'ci/jenkins/jenkinsfile.stellar-ai-cpu')
job.definition.lightweight = true
def retention = new LogRotator(-1, 2, -1, 2)
retention.removeLastBuild = true
job.buildDiscarder = retention
job.addProperty(new DisableConcurrentBuildsJobProperty())
job.addProperty(new ParametersDefinitionProperty([
    new StringParameterDefinition('CANDIDATE_PR', '', 'GitHub PR number; otherwise supply candidate and baseline refs.'),
    new StringParameterDefinition('CANDIDATE_REF', '', 'Candidate branch or full commit SHA.'),
    new StringParameterDefinition('BASELINE_REF', '', 'Baseline branch or full commit SHA.'),
    new BooleanParameterDefinition('ALLOW_BEHIND_CANDIDATE', false, 'Allow historical comparisons.')
]))
job.save()
j.setInstallState(InstallState.INITIAL_SETUP_COMPLETED)
j.save()
marker.text = 'Configured; this hook will not overwrite later administrator changes.\n'
println('Stellar AI CPU configured. Admin password and CLI token are in private files in JENKINS_HOME.')
