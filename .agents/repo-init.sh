#!/bin/bash

# Determine project root relative to the script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Define paths relative to project root
SRC_DIR="$PROJECT_ROOT/.agents/skills"
CLAUDE_DEST="$PROJECT_ROOT/.claude/skills"
CODEX_DEST="$PROJECT_ROOT/.codex/skills"

echo "🔄 Syncing Agent Skills across Claude and Codex environments..."

# Ensure target directories exist
mkdir -p "$CLAUDE_DEST"
mkdir -p "$CODEX_DEST"

# Each skill lives in its own directory: .agents/skills/<name>/SKILL.md
# Mirror each skill directory (as a whole) into .claude/skills/<name> and
# .codex/skills/<name>, since both tools discover skills at that path shape.
if [ -d "$SRC_DIR" ]; then
    for skill_dir in "$SRC_DIR"/*/; do
        # Ensure it's a directory
        [ -d "$skill_dir" ] || continue

        skill_name=$(basename "$skill_dir")

        # Create relative symlinks to prevent broken paths on different machines.
        # -n keeps this idempotent: it replaces an existing symlink instead of
        # nesting a new one inside an already-linked directory.
        ln -sfn "../../.agents/skills/$skill_name" "$CLAUDE_DEST/$skill_name"
        ln -sfn "../../.agents/skills/$skill_name" "$CODEX_DEST/$skill_name"
    done
    echo "✅ Symlinks successfully generated in .claude/skills/ and .codex/skills/"
else
    echo "❌ Error: Source directory $SRC_DIR does not exist."
    exit 1
fi

