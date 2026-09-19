#!/usr/bin/env bash
# Run after authenticating the computer to the intended tailnet.
set -euo pipefail
(( EUID == 0 )) || { echo 'Run with sudo to configure Cockpit origins.' >&2; exit 1; }
tailnet_host=$(tailscale status --json | python3 -c 'import json,sys; s=json.load(sys.stdin); assert s["BackendState"] == "Running", "Sign in with tailscale up first"; print(s["Self"]["DNSName"].rstrip("."))')
[[ "$tailnet_host" =~ ^[a-zA-Z0-9.-]+\.ts\.net$ ]] || { echo 'No valid tailnet DNS name found.' >&2; exit 1; }
python3 - "$tailnet_host" <<'PY'
from configparser import ConfigParser
from pathlib import Path
import sys
path = Path('/etc/cockpit/cockpit.conf')
config = ConfigParser(interpolation=None)
config.optionxform = str
config.read(path)
config['WebService']['Origins'] = 'http://localhost:9090 http://127.0.0.1:9090 https://' + sys.argv[1] + ':8443'
with path.open('w') as stream:
    config.write(stream)
PY
systemctl restart cockpit-container.service
tailscale serve --bg --https=443 http://127.0.0.1:8080
tailscale serve --bg --https=8443 http://127.0.0.1:9090
printf 'Jenkins: https://%s/\nCockpit: https://%s:8443/\n' "$tailnet_host" "$tailnet_host"
