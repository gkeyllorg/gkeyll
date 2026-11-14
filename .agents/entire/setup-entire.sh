#!/usr/bin/env bash
#
# First-time Entire setup for Gkeyll.
#
# Usage:
#   ./.agents/entire/setup-entire.sh [options]
#
# Prerequisite:
#   Install the Entire CLI before running this script:
#   https://docs.entire.io/installation
#
# Options:
#   --repo REPO       Repository used to store Entire checkpoint data.
#                     Default: gkeyll-agvcs
#
#   --agent AGENTS    Agent(s) to configure. Multiple agents may be
#                     comma-separated, e.g. codex,claude-code.
#                     If omitted, Entire prompts interactively.
#
#   -h, --help        Show this help.
#
# Examples:
#   ./.agents/entire/setup-entire.sh
#   ./.agents/entire/setup-entire.sh --agent codex
#   ./.agents/entire/setup-entire.sh --agent codex,claude-code
#   ./.agents/entire/setup-entire.sh --repo gkeyll --agent codex
#

set -euo pipefail

ORG="gkeyllorg"
REPO="gkeyll-agvcs"
AGENTS=""

usage()
{
  cat <<EOF
Usage:
  $0 [options]

Prerequisite:
  Install the Entire CLI before running this script:
  https://docs.entire.io/installation

Options:
  --repo REPO       Repository used to store Entire checkpoint data.
                    Default: gkeyll-agvcs

  --agent AGENTS    Agent(s) to configure. Multiple agents may be
                    comma-separated, e.g. codex,claude-code.
                    If omitted, Entire prompts interactively.

  -h, --help        Show this help.

Examples:
  $0
  $0 --agent codex
  $0 --agent codex,claude-code
  $0 --repo gkeyll --agent codex
EOF
}

fail()
{
  echo "ERROR: $*" >&2
  exit 1
}

# ----------------------------------------------------------------------
# Parse arguments
# ----------------------------------------------------------------------

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo)
      [[ $# -ge 2 ]] || fail "--repo requires an argument."
      REPO="$2"
      shift 2
      ;;
    --agent)
      [[ $# -ge 2 ]] || fail "--agent requires an argument."
      AGENTS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      fail "Unknown option: $1. Use --help for usage."
      ;;
  esac
done

[[ "$REPO" == */* ]] || REPO="$ORG/$REPO"

# ----------------------------------------------------------------------
# Basic checks
# ----------------------------------------------------------------------

command -v git >/dev/null 2>&1 \
  || fail "git is not installed."

command -v entire >/dev/null 2>&1 \
  || fail "Entire is not installed or is not in PATH."

ROOT="$(git rev-parse --show-toplevel 2>/dev/null)" \
  || fail "Run this script from inside a Gkeyll checkout."

cd "$ROOT"

# This script is intended for first-time setup.
if [[ -e .entire/settings.local.json ]] ||
   git for-each-ref --format='%(refname)' \
     refs/entire/checkpoints/ | grep -q . ||
   git show-ref --verify --quiet \
     refs/heads/entire/checkpoints/v1; then

  cat >&2 <<'EOF'
ERROR: Existing Entire configuration or checkpoint data was found.

If you want to change where existing Entire data is stored, use:

  ./.agents/entire/migrate-entire.sh

No changes were made.
EOF
  exit 1
fi

# ----------------------------------------------------------------------
# Construct a URL for another GitHub repository using the same
# authentication/transport style as origin.
#
# Examples supported:
#
#   https://github.com/user/gkeyll.git
#   https://janedoe@github.com/user/gkeyll.git
#   git@github.com:user/gkeyll.git
#   git@github-work:user/gkeyll.git
#   ssh://git@github.com/user/gkeyll.git
# ----------------------------------------------------------------------

ORIGIN="$(git remote get-url origin 2>/dev/null || true)"
[[ -n "$ORIGIN" ]] || fail "No Git remote named 'origin' was found."

case "$ORIGIN" in
  git@*:*)
    SSH_HOST="${ORIGIN%%:*}"
    URL="${SSH_HOST}:${REPO}.git"
    ;;

  ssh://*)
    REST="${ORIGIN#ssh://}"
    SSH_AUTHORITY="${REST%%/*}"
    URL="ssh://${SSH_AUTHORITY}/${REPO}.git"
    ;;

  http://*github.com/*|https://*github.com/*)
    if [[ "$ORIGIN" =~ ^(https?://)([^/]+@)?github\.com/ ]]; then
      HTTPS_BASE="${BASH_REMATCH[1]}${BASH_REMATCH[2]:-}github.com"
      URL="${HTTPS_BASE}/${REPO}.git"
    else
      fail "Could not parse GitHub HTTPS origin: $ORIGIN"
    fi
    ;;

  *)
    fail "Could not determine how to access GitHub from origin: $ORIGIN"
    ;;
esac

# ----------------------------------------------------------------------
# Make public storage an explicit choice.
# ----------------------------------------------------------------------

if [[ "$REPO" == "$ORG/gkeyll" ]]; then
  echo
  echo "WARNING: Entire checkpoint/session data stored in"
  echo
  echo "  $REPO"
  echo
  echo "will be PUBLIC."
  echo

  read -r -p "Continue? [y/N] " reply
  [[ "$reply" =~ ^[Yy]$ ]] || exit 1
fi

# ----------------------------------------------------------------------
# Verify read/write access.
# ----------------------------------------------------------------------

echo "Checking access to $REPO..."

if ! git ls-remote "$URL" >/dev/null 2>&1; then
  echo
  echo "ERROR: Cannot access $REPO." >&2

  if [[ "$REPO" == "$ORG/gkeyll-agvcs" ]]; then
    cat >&2 <<EOF

gkeyll-agvcs is private. Ask a gkeyllorg administrator for
read/write access if you would like to use it.

You may still develop Gkeyll normally without Entire.

If you intentionally want to store Entire data in the PUBLIC
Gkeyll repository for now, run:

  ./.agents/entire/setup-entire.sh --repo gkeyll
EOF
  fi

  exit 1
fi

# Check write permission without actually creating the ref.
git push --dry-run --no-verify \
  "$URL" HEAD:refs/entire/access-check/setup \
  >/dev/null 2>&1 \
  || fail "Write access to $REPO is required."

# ----------------------------------------------------------------------
# Enable Entire
# ----------------------------------------------------------------------

echo
echo "Configuring Entire checkpoint storage:"
echo
echo "  $REPO"
echo

if [[ -z "$AGENTS" ]]; then

  # Let Entire prompt interactively for agent selection.
  entire enable \
    --local \
    --checkpoint-remote "github:$REPO"

else

  IFS=',' read -r -a agent_list <<< "$AGENTS"
  cleaned_agents=()

  for agent in "${agent_list[@]}"; do
    agent="${agent//[[:space:]]/}"
    [[ -n "$agent" ]] \
      || fail "Empty agent name in --agent '$AGENTS'."
    cleaned_agents+=("$agent")
  done

  first_agent="${cleaned_agents[0]}"

  echo "Configuring agent: $first_agent"

  entire enable \
    -y \
    --local \
    --agent "$first_agent" \
    --checkpoint-remote "github:$REPO"

  # Add any additional agents after initial setup.
  if [[ ${#cleaned_agents[@]} -gt 1 ]]; then
    for agent in "${cleaned_agents[@]:1}"; do
      echo "Adding agent: $agent"
      entire agent add "$agent"
    done
  fi
fi

echo
echo "Entire setup complete."
echo
entire status
