#!/usr/bin/env bash
# Apply the host settings needed by the user services. Run with sudo.
set -euo pipefail
(( EUID == 0 )) || { echo 'Run with sudo; these are system service settings.' >&2; exit 1; }
target_user=${1:-${SUDO_USER:-}}
[[ -n "$target_user" && "$target_user" != root ]] || { echo 'Specify the desktop username.' >&2; exit 1; }
id "$target_user" >/dev/null

# Suspend stops the CPU, GPU and Jenkins. Display blanking remains available.
systemctl mask sleep.target suspend.target hibernate.target hybrid-sleep.target suspend-then-hibernate.target
loginctl enable-linger "$target_user"
systemctl enable --now tailscaled.service
tailscale set --operator="$target_user"

# Cockpit's current container authenticates to the host over local SSH.
# A loopback listener supplies its terminal without opening a network SSH port.
install -d -m 755 /etc/ssh/sshd_config.d
cat > /etc/ssh/sshd_config.d/00-gkeyll-cockpit.conf <<'EOF'
ListenAddress 127.0.0.1
ListenAddress ::1
Match Address 127.0.0.1,::1
    PasswordAuthentication yes
Match all
EOF
ssh-keygen -A
/usr/sbin/sshd -t
systemctl enable --now sshd.service
systemctl reload sshd.service

# Use Bazzite's existing Cockpit container, bound only to loopback.
install -d -m 755 /etc/containers/systemd/cockpit-container.container.d
cat > /etc/containers/systemd/cockpit-container.container.d/20-gkeyll-loopback.conf <<'EOF'
[Container]
Image=quay.io/cockpit/ws@sha256:f243418febe14ea7bb4f7d76ac4a4c72093d81ac4bf68b2270670b100f8cc323
Exec=/container/label-run --address=127.0.0.1 --no-tls

[Service]
TimeoutStartSec=600
EOF
install -d -m 755 /etc/cockpit
python3 - <<'PY'
from configparser import ConfigParser
from pathlib import Path
path = Path('/etc/cockpit/cockpit.conf')
config = ConfigParser(interpolation=None)
config.optionxform = str
if path.exists():
    backup = path.with_suffix('.conf.before-gkeyll')
    if not backup.exists():
        backup.write_bytes(path.read_bytes())
    config.read(path)
if not config.has_section('WebService'):
    config.add_section('WebService')
config['WebService']['AllowUnencrypted'] = 'true'
config['WebService']['ProtocolHeader'] = 'X-Forwarded-Proto'
config['WebService']['Origins'] = 'http://localhost:9090 http://127.0.0.1:9090'
with path.open('w') as stream:
    config.write(stream)
PY
systemctl daemon-reload
systemctl enable --now cockpit.service
systemctl restart cockpit-container.service
echo 'Host ready: suspend disabled, lingering enabled, Tailscale and loopback Cockpit started.'
echo 'Next: sign in with tailscale up, then run setup-tailnet.sh with sudo.'
