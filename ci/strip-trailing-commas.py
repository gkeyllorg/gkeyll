#!/usr/bin/env python3
"""Strip trailing commas immediately before a closing '}' in C/C++/CUDA source.

clang-format expands every element of a braced initializer list onto its own
line whenever the list already ends with a trailing comma, regardless of
column limit. Since some authors habitually add one and others don't, the
same kind of initializer formats inconsistently depending on who last edited
it. This script removes the comma before it ever reaches clang-format, so
formatting stays deterministic.

Usage:
    strip-trailing-commas.py FILE...          # rewrite files in place
    strip-trailing-commas.py --check FILE...  # report only; exit 1 if any file would change
"""
import re
import sys

# Matches an (optional) encoding prefix + R"delim(...)delim" raw string literal.
# The delimiter is captured and reused via a backreference so the literal's
# body -- which may itself contain quotes, braces, or comment-like sequences
# -- is never misread as comment/string/brace syntax.
_RAW_STRING_RE = re.compile(
    r'(?:u8|u|U|L)?R"([^\s"\\()]{0,16})\((?:(?!\)\1").)*\)\1"', re.DOTALL
)

_TRAILING_COMMA_RE = re.compile(r",\s*\}")


def _blank(out, start, end, filler):
    for k in range(start, end):
        if out[k] != "\n":
            out[k] = filler


def _mask(text):
    """Return a same-length copy of text safe to run the trailing-comma regex on.

    Comment contents are blanked to spaces, since a comment between a comma
    and a '}' is trivia -- \\s* should swallow it. String/char/raw-string
    literal contents are blanked to 'x' instead: they're real code (e.g. an
    array element that happens to be a string literal), so they must NOT
    look like whitespace, or a comma several elements back could wrongly
    appear to be immediately followed by '}'. Delimiters and newlines are
    preserved either way so offsets line up with the original text.
    """
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
    comma_offsets = [m.start() for m in _TRAILING_COMMA_RE.finditer(masked)]
    if not comma_offsets:
        return text
    chars = list(text)
    for offset in reversed(comma_offsets):
        del chars[offset]
    return "".join(chars)


def main(argv):
    check = "--check" in argv
    files = [a for a in argv if a != "--check"]
    changed = []
    for path in files:
        with open(path, "r") as f:
            original = f.read()
        updated = strip_trailing_commas(original)
        if updated != original:
            changed.append(path)
            if not check:
                with open(path, "w") as f:
                    f.write(updated)
    if check:
        for path in changed:
            print(f"trailing comma(s) found: {path}")
        return 1 if changed else 0
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
