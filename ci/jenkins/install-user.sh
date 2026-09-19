#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
(( EUID != 0 )) || { echo 'Run this as your desktop user, without sudo.' >&2; exit 1; }
state="$HOME/.local/share/gkeyll-jenkins"
units="$HOME/.config/containers/systemd"
umask 077
mkdir -p "$state" "$state/controller" "$state/agent" "$state/agent-credentials" "$units" "$HOME/.config/gkeyll"
if [[ ! -s "$state/admin-password" ]]; then
  python3 -c 'import secrets; print(secrets.token_urlsafe(32))' > "$state/admin-password"
fi
podman build -f Containerfile.controller -t localhost/gkeyll-jenkins:local .
podman build -f Containerfile.gpu -t localhost/gkeyll-gpu-agent:local .
# Jenkins preserves initialization files in its data directory across upgrades.
install -D -m 644 bootstrap.groovy "$state/controller/init.groovy.d/gkeyll.groovy"
install -m 644 gkeyll-jenkins.network gkeyll-jenkins.container gkeyll-gpu-agent.container "$units/"
systemctl --user daemon-reload
systemctl --user restart gkeyll-jenkins.service gkeyll-gpu-agent.service
if [[ ! -x "$state/java/bin/java" ]]; then
  podman cp gkeyll-jenkins:/opt/java/openjdk "$state/java"
fi
cat > "$HOME/.config/gkeyll/jenkins.env" <<'EOF'
export JENKINS_URL=http://127.0.0.1:8080
export JENKINS_CLI_AUTH_FILE="$HOME/.config/gkeyll/jenkins-cli.auth"
export JAVA_HOME="$HOME/.local/share/gkeyll-jenkins/java"
export PATH="$JAVA_HOME/bin:$PATH"
EOF
printf 'Jenkins: http://localhost:8080/\nUsername: mrosen\nPassword file: %s/admin-password\n' "$state"
