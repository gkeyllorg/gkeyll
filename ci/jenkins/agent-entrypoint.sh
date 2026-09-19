#!/usr/bin/env bash
set -euo pipefail
umask 077
: "${JENKINS_URL:=http://gkeyll-jenkins:8080/}"
: "${JENKINS_AGENT_NAME:=bazzite-rtx3090}"
until [[ -s /run/agent-credentials/secret ]] && \
  curl --fail --silent --show-error --connect-timeout 5 --max-time 30 \
    "${JENKINS_URL%/}/jnlpJars/agent.jar" -o agent.jar.tmp; do
  sleep 10
done
mv agent.jar.tmp agent.jar
exec java -Xmx512m -jar agent.jar -url "$JENKINS_URL" \
  -name "$JENKINS_AGENT_NAME" -secret @/run/agent-credentials/secret \
  -webSocket -workDir /home/jenkins/agent
