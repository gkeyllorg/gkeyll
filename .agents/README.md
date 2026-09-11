# Agent context

Maintain project-wide instructions in [`AGENTS.md`](../AGENTS.md) and detailed
procedures in `skills/<skill-name>/SKILL.md`. The root instructions route each
task to the relevant procedure, so agents need not load every skill at startup.

## Tool entry points

| Tool | Entry point | Detailed context |
| --- | --- | --- |
| Codex | `AGENTS.md` | Native discovery of `.agents/skills/` |
| Claude Code | `CLAUDE.md` imports `AGENTS.md` | Read skills through the shared task table |
| Gemini CLI | `GEMINI.md` imports `AGENTS.md` | Read skills through the shared task table |
| GitHub Copilot | `.github/copilot-instructions.md` points to `AGENTS.md` | Read skills through the shared task table in modes with repository file access |
| Other agents | Configure the tool to read `AGENTS.md` | Read the referenced Markdown files |

The entry points are regular text files. No symlink support, setup script,
generated copies, or user-specific absolute paths are required. Claude's native
project slash commands are not registered by this layout; the same procedures
are available through task routing or an explicit request to read a skill.
Copilot's pointer is an instruction to read a file, not a native import; support
depends on the Copilot surface and its access to repository files.

Keep tool-specific configuration in separate real directories (`.claude/`,
`.codex/`, etc.) if needed. Do not alias those directories to `.agents/`:
their settings have different meanings, and protected-directory symlinks can
prevent sandbox startup. Local settings and permissions belong in ignored files
or the tool's user configuration, not in shared skills.

## Maintaining context

- Edit shared policy only in `AGENTS.md`; keep tool entry points minimal.
- Use lowercase hyphenated skill names, matching the containing directory.
  Each skill needs YAML frontmatter with `name` and a `description` that says
  when it applies. Keep procedures in the body.
- When adding or renaming a skill, update the task table and relative links.
- Keep machine-specific build paths in the existing machine configuration;
  skills should explain how to find the appropriate configuration.
- Add an adapter for another tool only after checking its documented discovery
  mechanism. `AGENTS.md` and `@` imports are not universal conventions.

## Verification

After changing this layout, check that every local Markdown link resolves, skill
names match their folders, frontmatter is valid, and entry points are regular
files. Run `git diff --check`. A fresh checkout should need no setup step.

Start a fresh session in each tool you use and ask:

> Which repository instructions apply, and which skill files would you read
> before changing a gyrokinetic updater and running its tests? Do not edit files
> or run simulations.

Expect the shared policy and relevant design, naming, gyrokinetic, compilation,
and testing guidance. Also check an unrelated documentation task to ensure it
does not load every solver skill. Claude's `/context` and Gemini's `/memory show`
can help inspect loaded project instructions. Static link checks alone do not
prove that a particular client version follows the routing.

## Discovery references

- [Codex instructions](https://developers.openai.com/codex/guides/agents-md/)
  and [skills](https://developers.openai.com/codex/skills/)
- [Claude Code instructions and imports](https://code.claude.com/docs/en/memory)
- [Gemini CLI context and imports](https://geminicli.com/docs/cli/gemini-md/)
- [Copilot repository instructions](https://docs.github.com/en/copilot/how-tos/configure-custom-instructions-in-your-ide/add-repository-instructions-in-your-ide)
