#!/usr/bin/env python3
"""Format Gkeyll C/C++/CUDA sources: clang-format plus the trailing-comma rules.

Each file goes through four steps, in memory, and is written back only if the
final result differs from what is on disk (so the hook is idempotent and a
clean tree passes):

  1. strip every trailing comma before a closing '}'   (author-independent input)
  2. clang-format
  3. add a trailing comma to every designated-initializer list that clang-format
     could not fit on one line, and to any braced list it left starting a line
     with '){' (a compound literal whose cast parentheses got broken instead)
  4. clang-format again -- the trailing comma is what makes it put one member
     per line, with the first member on its own line below the '{'

The net rule is: a designated initializer that fits on one line stays on one
line; one that does not gets one member per line. That is independent of
whether the author typed a trailing comma.

Usage:
    format.py FILE...          # rewrite files in place
    format.py --check FILE...  # report only; exit 1 if any file would change

Requires clang-format (the pinned version, see .pre-commit-config.yaml) on PATH.
"""
import re
import subprocess
import sys

# Matches an (optional) encoding prefix + R"delim(...)delim" raw string literal.
_RAW_STRING_RE = re.compile(
    r'(?:u8|u|U|L)?R"([^\s"\\()]{0,16})\((?:(?!\)\1").)*\)\1"', re.DOTALL
)
_TRAILING_COMMA_RE = re.compile(r",\s*\}")
# First token of a designated initializer: '.name' or '[idx]' ('.5' is a float).
_DESIGNATOR_RE = re.compile(r"\.[A-Za-z_]|\[")


def _blank(out, start, end, filler):
    for k in range(start, end):
        if out[k] != "\n":
            out[k] = filler


def _mask(text):
    """Same-length copy of text with comments blanked to spaces and string,
    char and raw-string literal contents blanked to 'x', so brace/comma
    scanning never sees syntax inside them. Newlines are kept."""
    out = list(text)
    i = 0
    n = len(text)
    while i < n:
        two = text[i : i + 2]
        if two == "//":
            j = text.find("\n", i)
            end = n if j == -1 else j
            _blank(out, i, end, " ")
            i = end
        elif two == "/*":
            end = text.find("*/", i + 2)
            end = n if end == -1 else end + 2
            _blank(out, i, end, " ")
            i = end
        elif text[i] in "'\"":
            quote = text[i]
            j = i + 1
            while j < n and text[j] != quote:
                if text[j] == "\\":
                    j += 1
                j += 1
            end = min(j + 1, n)
            _blank(out, i, end, "x")
            i = end
        elif text[i] in "uUL8R" and _RAW_STRING_RE.match(text, i):
            end = _RAW_STRING_RE.match(text, i).end()
            _blank(out, i, end, "x")
            i = end
        else:
            i += 1
    return "".join(out)


def strip_trailing_commas(text):
    masked = _mask(text)
    offsets = [m.start() for m in _TRAILING_COMMA_RE.finditer(masked)]
    if not offsets:
        return text
    chars = list(text)
    for off in reversed(offsets):
        del chars[off]
    return "".join(chars)


def add_trailing_commas(text):
    m = _mask(text)
    inserts = []
    stack = []
    for i, ch in enumerate(m):
        if ch == "{":
            stack.append(i)
        elif ch == "}" and stack:
            o = stack.pop()
            body = m[o + 1 : i]
            designated = bool(_DESIGNATOR_RE.match(body.lstrip()))
            line_start = m.rfind("\n", 0, o) + 1
            cast_broken = m[line_start:o].strip() == ")"
            if not (cast_broken or (designated and "\n" in body)):
                continue
            j = i - 1
            # Skip whitespace and the backslash of a macro line continuation.
            while j > o and (m[j].isspace() or m[j] == "\\"):
                j -= 1
            if m[j] != ",":
                inserts.append(j + 1)
    if not inserts:
        return text
    chars = list(text)
    for off in sorted(inserts, reverse=True):
        chars.insert(off, ",")
    return "".join(chars)


def clang_format(text, path):
    res = subprocess.run(
        ["clang-format", "--style=file", f"--assume-filename={path}"],
        input=text, capture_output=True, text=True, check=True,
    )
    return res.stdout


def format_text(text, path):
    text = strip_trailing_commas(text)
    text = clang_format(text, path)
    text = add_trailing_commas(text)
    return clang_format(text, path)


def main(argv):
    check = "--check" in argv
    files = [a for a in argv if a != "--check"]
    changed = []
    for path in files:
        with open(path, "r") as f:
            original = f.read()
        updated = format_text(original, path)
        if updated != original:
            changed.append(path)
            if not check:
                with open(path, "w") as f:
                    f.write(updated)
    if check:
        for path in changed:
            print(f"not formatted: {path}")
        return 1 if changed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
