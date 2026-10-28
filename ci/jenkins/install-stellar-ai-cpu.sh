#!/usr/bin/env bash
# Install a new, account-owned controller in shared storage. No sudo required.
set -euo pipefail
umask 077
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/../.." && pwd)"
export GKEYLL_CI_ROOT="${GKEYLL_CI_ROOT:-/scratch/gpfs/CSD/${USER:?}/gkeyll_ci_stellar_ai_cpu}"
export JAVA_HOME="${JAVA_HOME:-/usr/lib/jvm/jre-21-openjdk}"
root="$GKEYLL_CI_ROOT"
[[ "$root" = /* && "$root" != / && "$root" != *..* ]] || { echo 'Use an absolute shared scratch directory' >&2; exit 1; }
[[ ! -e "$root/jenkins_home/config.xml" ]] || { echo 'Controller already configured; refusing to overwrite it' >&2; exit 1; }
[[ -x "$JAVA_HOME/bin/java" ]] || { echo 'Java 21 is required' >&2; exit 1; }
command -v tmux >/dev/null
mkdir -p "$root/trusted/ci/jenkins" "$root/trusted/machines" "$root/jenkins_home/init.groovy.d"
chmod 700 "$root" "$root/jenkins_home"
curl -fLsS --retry 2 -o "$root/jenkins.war.tmp" https://get.jenkins.io/war-stable/2.580.1/jenkins.war
printf '%s  %s\n' 393bf2476352dd726519fd1f92ce66eac7d23c0936e6967420b717060c4c40d0 "$root/jenkins.war.tmp" | sha256sum --check
mv "$root/jenkins.war.tmp" "$root/jenkins.war"
curl -fLsS --retry 2 -o "$root/plugin-manager.jar" https://github.com/jenkinsci/plugin-installation-manager-tool/releases/download/2.15.0/jenkins-plugin-manager-2.15.0.jar
"$JAVA_HOME/bin/java" -jar "$root/plugin-manager.jar" --war "$root/jenkins.war" \
  --plugin-download-directory "$root/jenkins_home/plugins" \
  --plugins workflow-aggregator git credentials-binding
for name in jenkins-stellar-ai-cpu.sh jenkins-stellar_cpu.sh jenkins-client-artifacts.sh gkeyll-ci.sh; do
  cp "$script_dir/$name" "$root/trusted/ci/jenkins/"
done
revision="$(git -C "$repo_root" rev-parse HEAD)"
digest="$(cd "$root/trusted" && find ci -type f ! -name deployment-revision.txt -print0 | sort -z | xargs -0 sha256sum | sha256sum | cut -d' ' -f1)"
printf '%s+snapshot.%s\n' "$revision" "$digest" > "$root/trusted/ci/jenkins/deployment-revision.txt"
cp "$script_dir/setup-stellar-ai-cpu.groovy" "$root/jenkins_home/init.groovy.d/10-stellar-ai-cpu.groovy"
cp "$script_dir/queue-status.groovy" "$root/jenkins_home/init.groovy.d/20-gkeyll-queue-status.groovy"
"$root/trusted/ci/jenkins/gkeyll-ci.sh" stellar-ai-cpu start
printf 'Login: admin; password file: %s/jenkins_home/admin-password\n' "$root"
