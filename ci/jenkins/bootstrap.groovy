import jenkins.model.Jenkins
import jenkins.model.JenkinsLocationConfiguration
import jenkins.install.InstallState
import jenkins.slaves.JnlpAgentReceiver
import hudson.model.Node
import hudson.security.HudsonPrivateSecurityRealm
import hudson.security.FullControlOnceLoggedInAuthorizationStrategy
import hudson.security.GlobalMatrixAuthorizationStrategy
import hudson.model.Item
import org.jenkinsci.plugins.matrixauth.PermissionEntry
import org.jenkinsci.plugins.matrixauth.AuthorizationType
import hudson.security.csrf.DefaultCrumbIssuer
import hudson.slaves.DumbSlave
import hudson.slaves.JNLPLauncher
import hudson.slaves.RetentionStrategy
import org.jenkinsci.plugins.workflow.job.WorkflowJob
import org.jenkinsci.plugins.workflow.cps.CpsFlowDefinition
import hudson.slaves.EnvironmentVariablesNodeProperty
import hudson.model.ParametersDefinitionProperty
import hudson.model.StringParameterDefinition
import hudson.model.User
import jenkins.security.ApiTokenProperty
import java.nio.file.Files
import java.nio.file.attribute.PosixFilePermissions
import org.jenkinsci.plugins.workflow.multibranch.WorkflowMultiBranchProject
import org.jenkinsci.plugins.pipeline.multibranch.defaults.PipelineBranchDefaultsProjectFactory
import org.jenkinsci.plugins.configfiles.GlobalConfigFiles
import org.jenkinsci.plugins.configfiles.groovy.GroovyScript
import org.jenkinsci.plugins.github_branch_source.GitHubSCMSource
import org.jenkinsci.plugins.github_branch_source.BranchDiscoveryTrait
import org.jenkinsci.plugins.github_branch_source.OriginPullRequestDiscoveryTrait
import org.jenkinsci.plugins.github_branch_source.ForkPullRequestDiscoveryTrait
import jenkins.scm.impl.trait.WildcardSCMHeadFilterTrait
import jenkins.branch.BranchSource
import com.cloudbees.hudson.plugins.folder.computed.PeriodicFolderTrigger
import com.cloudbees.hudson.plugins.folder.computed.DefaultOrphanedItemStrategy

def j = Jenkins.get()
def marker = new File(j.rootDir, '.gkeyll-initialized')
if (!marker.exists()) {
    def realm = new HudsonPrivateSecurityRealm(false)
    realm.createAccount('mrosen', new File('/run/secrets/admin-password').text.trim())
    j.setSecurityRealm(realm)
    def auth = new FullControlOnceLoggedInAuthorizationStrategy()
    auth.setAllowAnonymousRead(false)
    j.setAuthorizationStrategy(auth)
    j.setCrumbIssuer(new DefaultCrumbIssuer(true))
    j.setInstallState(InstallState.INITIAL_SETUP_COMPLETED)
    j.setNumExecutors(0)
    j.setSlaveAgentPort(-1)
    JenkinsLocationConfiguration.get().setUrl('http://localhost:8080/')
    j.save()
    marker.text = 'Initial security configuration completed.\n'
}
// Keep administration with the owner; accounts created for colleagues can
// inspect, submit, and cancel builds without changing scripts or credentials.
if (j.authorizationStrategy instanceof FullControlOnceLoggedInAuthorizationStrategy) {
    def auth = new GlobalMatrixAuthorizationStrategy()
    auth.add(Jenkins.ADMINISTER, new PermissionEntry(AuthorizationType.USER, 'mrosen'))
    def members = new PermissionEntry(AuthorizationType.GROUP, 'authenticated')
    [Jenkins.READ, Item.READ, Item.BUILD, Item.CANCEL].each { auth.add(it, members) }
    j.setAuthorizationStrategy(auth)
}
def name = 'bazzite-rtx3090'
if (j.getNode(name) == null) {
    def node = new DumbSlave(name, '/home/jenkins/agent', new JNLPLauncher())
    node.setNumExecutors(1)
    node.setLabelString('linux bazzite cuda rtx3090')
    node.setMode(Node.Mode.EXCLUSIVE)
    node.setRetentionStrategy(new RetentionStrategy.Always())
    j.addNode(node)
}
def secret = new File('/run/agent-credentials/secret')
secret.text = JnlpAgentReceiver.DATABASE.getSecretOf(name)
Files.setPosixFilePermissions(secret.toPath(), PosixFilePermissions.fromString('rw-------'))
// Install a trusted snapshot of the repository's existing personal Pipeline.
// Candidate PRs supply simulation code, never the controlling Pipeline.
def pipelineDir = new File('/usr/local/share/gkeyll')
def personal = new File(pipelineDir, 'jenkinsfile.personal').text.replace(
    "readTrusted('ci/jenkins/check_regression_results.lua')",
    new File(pipelineDir, 'check_regression_results.lua').text.inspect())
def job = j.getItem('gkeyll-ci-personal')
if (job == null) {
    job = j.createProject(WorkflowJob, 'gkeyll-ci-personal')
    job.addProperty(new ParametersDefinitionProperty(
        new StringParameterDefinition('CANDIDATE_PR', '', 'PR number, or leave empty and set both references.'),
        new StringParameterDefinition('CANDIDATE_REF', '', 'Candidate branch or full commit SHA.'),
        new StringParameterDefinition('BASELINE_REF', '', 'Baseline branch or full commit SHA.')))
}
job.setDescription('Repository personal CI on the Bazzite RTX 3090. Select a PR or candidate/baseline comparison. Uses the trusted jenkinsfile.personal installed with the controller image.')
job.setDefinition(new CpsFlowDefinition(personal, true))
job.save()
def envProperty = j.globalNodeProperties.get(EnvironmentVariablesNodeProperty)
if (envProperty == null) {
    envProperty = new EnvironmentVariablesNodeProperty()
    j.globalNodeProperties.add(envProperty)
}
def defaults = [
    PERSONAL_NODE_LABEL: 'bazzite && cuda && rtx3090',
    PERSONAL_MKDEPS_SCRIPT: 'mkdeps.linux.sh',
    PERSONAL_CONFIGURE_SCRIPT: 'configure.linux.gpu.sh',
    PERSONAL_GITHUB_CREDENTIAL_ID: 'gkeyll-github',
    PERSONAL_BUILD_JOBS: '8', PERSONAL_REGRESSION_JOBS: '1',
    TEAM_WORKSTATION_NODE_LABEL: 'bazzite && cuda && rtx3090',
    TEAM_WORKSTATION_MKDEPS_SCRIPT: 'mkdeps.linux.sh',
    TEAM_WORKSTATION_CONFIGURE_SCRIPT: 'configure.linux.gpu.sh',
    TEAM_WORKSTATION_GITHUB_CREDENTIAL_ID: 'gkeyll-github',
    TEAM_WORKSTATION_BUILD_JOBS: '8', TEAM_WORKSTATION_REGRESSION_JOBS: '1',
    MAKEFLAGS: 'CUDA_ARCH=86', OPENBLAS_NUM_THREADS: '1', OMP_NUM_THREADS: '1'
]
defaults.each { key, value -> envProperty.envVars.putIfAbsent(key, value) }
// A controller-managed default lets all existing PRs use the repository's
// team Pipeline, including branches which do not yet contain CI files.
def team = new File(pipelineDir, 'jenkinsfile.team_workstation').text.replace(
    "readTrusted('ci/jenkins/jenkinsfile.personal')", personal.inspect())
GlobalConfigFiles.get().save(new GroovyScript('gkeyll-team-pipeline',
    'Gkeyll team workstation', 'Trusted snapshot installed with the controller image.', team))
def multibranch = j.getItem('gkeyll-ci-team-workstation')
if (multibranch == null) {
    multibranch = j.createProject(WorkflowMultiBranchProject, 'gkeyll-ci-team-workstation')
    def source = new GitHubSCMSource('gkeyllorg', 'gkeyll')
    source.setId('gkeyll-public-source')
    // Authenticate API discovery; public Git checkouts in the Pipeline are
    // still anonymous. GitHub's anonymous quota cannot sustain team scans.
    source.setCredentialsId('gkeyll-github')
    source.setTraits([
        new BranchDiscoveryTrait(1),
        new OriginPullRequestDiscoveryTrait(2),
        new ForkPullRequestDiscoveryTrait(2, new ForkPullRequestDiscoveryTrait.TrustNobody()),
        new WildcardSCMHeadFilterTrait('main PR-*', '')
    ])
    multibranch.getSourcesList().add(new BranchSource(source))
    def factory = new PipelineBranchDefaultsProjectFactory()
    factory.setScriptId('gkeyll-team-pipeline')
    factory.setUseSandbox(true)
    multibranch.setProjectFactory(factory)
    multibranch.setOrphanedItemStrategy(new DefaultOrphanedItemStrategy(true, '14', '20'))
    multibranch.addTrigger(new PeriodicFolderTrigger('2m'))
    multibranch.setDescription('Shared Gkeyll RTX 3090 CI. Discovers main and PRs from all authors; the team Pipeline skips PRs targeting other branches. Uses the trusted controller Pipeline snapshot.')
    multibranch.save()
}
j.getDescriptorByType(hudson.tasks.Shell.DescriptorImpl).setShell('/bin/bash')
// The existing CLI expects a protected username:API-token file.
def cliAuth = new File('/run/cli-credentials/jenkins-cli.auth')
if (!cliAuth.exists()) {
    def user = User.getById('mrosen', false)
    def token = user.getProperty(ApiTokenProperty).tokenStore.generateNewToken('gkeyll-local-cli')
    user.save()
    cliAuth.text = 'mrosen:' + token.plainValue + '\n'
    Files.setPosixFilePermissions(cliAuth.toPath(), PosixFilePermissions.fromString('rw-------'))
}
j.save()
println('Gkeyll controller, GPU node, and team PR discovery configured.')
