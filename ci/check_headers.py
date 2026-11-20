#!/usr/bin/env python3
"""Check the header comment of regression inputs against the input itself.

Every C (creg/rt_*.c) and Lua (luareg/rt_*.lua) regression input starts with a
comment that describes the problem and states its figures of merit. A figure of
merit belongs to one resolution, one basis and one end time, so the header has
to state them, and they have to be those the input sets. For each input this
script checks:

  cells  The parenthesis after "Figures of merit" states cell counts, and they
         are the cell counts of the input (Nx, Ny, ..., then the velocity
         cells), all of them or the configuration-space ones.
  basis  A basis and polynomial order stated there, before any semicolon (what
         follows one describes a comparison run), are those of the input, and
         a Lua input named <test>_<basis>_p<order>.lua sets that basis and
         order.
  time   No time the header quotes as a plain number (t = 3.5, t_end = 8.5)
         lies beyond the end time of the input. Times in other units
         (t = 22/Omega_ci) are not checked.

An input whose header states no figure of merit is listed with --verbose and
is not an error.

Usage:
  ci/check_headers.py                 all inputs of vlasov (with the PKPM ones)
  ci/check_headers.py vlasov          all inputs of one solver
  ci/check_headers.py FILE...         the given inputs
  ci/check_headers.py --verbose ...   also list the inputs without a figure

The exit status is 1 if any header fails a check.
"""

import argparse
import glob
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SOLVERS = ("vlasov",)

# Names of the cell counts of an input, in the order a header lists them.
CELL_NAME = re.compile(r"^N(x|y|z|r|theta|phi|z1|z2|v\w*|p\w*)$")
C_ASSIGN = re.compile(r"^\s*(?:const\s+)?(?:int|double|bool)\s+(\w+)\s*=\s*([^;]+);", re.M)
LUA_ASSIGN = re.compile(r"^(?:local\s+)?(\w+)\s*=\s*([^\n]+?)\s*$", re.M)
FOM = re.compile(r"Figures? of merit[^()]*\(([^()]*)\)")
CELL_CLAIM = re.compile(r"(\d+(?:\s*x\s*\d+)+)\s*cells|(\d+)\s*cells")
BASIS_CLAIM = re.compile(r"(serendipity|tensor|hybrid)\s*p\s*=?\s*(\d)")
BASIS_SUFFIX = re.compile(r"_(ser|tensor|hyb)_p(\d)$")
PLAIN_TIME = re.compile(
  r"(?:\bto t\s*=|\bt_end\s*=|\bat t\s*=|\(t\s*=|\bover t\s*=|\bfor t\s*=)\s*"
  r"(\d+(?:\.\d*)?(?:e[-+]?\d+)?)(?![\d.]*\s*(?:/|pi\b|[A-Za-z_]))"
)
BASIS_NAME = {"ser": "serendipity", "tensor": "tensor", "hyb": "hybrid"}


def split_input(path):
  """The header (one string, comment markers removed) and the code of an input."""
  mark = "//" if path.endswith(".c") else "--"
  lines = open(path).read().split("\n")
  n = 0
  while n < len(lines) and lines[n].startswith(mark):
    n += 1
  header = " ".join(l[len(mark):].strip() for l in lines[:n])
  code = "\n".join(l.split(mark)[0] for l in lines[n:])
  return re.sub(r"\s+", " ", header), code


def scalar_values(path, code):
  """Scalars the input assigns, in order: name -> (text, number or None).

  A number is found for a literal and for an expression of earlier scalars.
  """
  rx = C_ASSIGN if path.endswith(".c") else LUA_ASSIGN
  env = {"pi": math.pi, "M_PI": math.pi, "sqrt": math.sqrt, "pow": pow, "exp": math.exp}
  values = {}
  for m in rx.finditer(code):
    name, text = m.group(1), m.group(2).strip().rstrip(",")
    if name in values:
      continue
    expr = text.replace("math.", "").replace("^", "**")
    try:
      number = float(eval(expr, {"__builtins__": {}}, env))  # noqa: S307 (inputs of this repository)
    except Exception:  # not a scalar expression
      number = None
    values[name] = (text, number)
    if number is not None:
      env[name] = number
  return values


def input_cells(values):
  """Cell counts of the input, in declaration order."""
  return [
    int(text) for name, (text, _) in values.items() if CELL_NAME.match(name) and text.isdigit()
  ]


def input_basis(path, values):
  """(basis, order) of the input; the basis is None when the input does not name one."""
  order = values.get("poly_order", (None, None))[1]
  basis = values.get("basis_type", (None, None))[0]
  if basis is not None:
    basis = basis.strip("\"'")
  return basis, None if order is None else int(order)


def check_input(path):
  """Failures of one input (strings), and whether its header states a figure of merit."""
  header, code = split_input(path)
  values = scalar_values(path, code)
  cells = input_cells(values)
  basis, order = input_basis(path, values)
  failures = []

  fom = FOM.search(header)
  has_figure = "igure of merit" in header or "igures of merit" in header
  if has_figure:
    claim = CELL_CLAIM.search(fom.group(1)) if fom else None
    if claim is None:
      failures.append("cells: the figures of merit state no cell count (input: %s)" % cells)
    else:
      stated = [int(x) for x in re.findall(r"\d+", claim.group(1) or claim.group(2))]
      if stated != cells and stated != cells[: len(stated)]:
        failures.append("cells: header %s, input %s" % (stated, cells))
    stated_basis = BASIS_CLAIM.search(fom.group(1).split(";")[0]) if fom else None
    if stated_basis:
      name, p = stated_basis.group(1), int(stated_basis.group(2))
      if order is not None and p != order:
        failures.append("basis: header p%d, input poly_order %d" % (p, order))
      if basis is not None and name != basis:
        failures.append("basis: header %s, input %s" % (name, basis))

  suffix = BASIS_SUFFIX.search(os.path.basename(path)[: -len(".lua")]) if path.endswith(".lua") else None
  if suffix:
    name, p = BASIS_NAME[suffix.group(1)], int(suffix.group(2))
    if order is not None and p != order:
      failures.append("basis: file name p%d, input poly_order %d" % (p, order))
    if basis is not None and name != basis:
      failures.append("basis: file name %s, input %s" % (name, basis))

  t_end = values.get("t_end", (None, None))[1]
  if t_end is not None:
    late = sorted({float(t) for t in PLAIN_TIME.findall(header) if float(t) > t_end * (1.0 + 1.0e-3)})
    if late:
      failures.append("time: header quotes t = %s, input ends at %g" % (", ".join("%g" % t for t in late), t_end))

  return failures, has_figure


def collect(args):
  """Input files named by the arguments (solvers, files; default: all solvers)."""
  files = []
  for a in args or SOLVERS:
    if a in SOLVERS:
      files += sorted(glob.glob(os.path.join(ROOT, a, "creg", "rt_*.c")))
      files += sorted(glob.glob(os.path.join(ROOT, a, "luareg", "rt_*.lua")))
    else:
      files.append(os.path.abspath(a))
  return files


def main():
  parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("inputs", nargs="*", help="solvers (vlasov) or input files")
  parser.add_argument("--verbose", action="store_true", help="list the inputs without a figure of merit")
  opts = parser.parse_args()

  files = collect(opts.inputs)
  num_failed, without_figure = 0, []
  for path in files:
    failures, has_figure = check_input(path)
    name = os.path.relpath(path, ROOT)
    if not has_figure:
      without_figure.append(name)
    for failure in failures:
      print("%s: %s" % (name, failure))
    num_failed += bool(failures)
  if opts.verbose:
    for name in without_figure:
      print("%s: no figure of merit" % name)
  print(
    "%d inputs checked, %d with a failure, %d without a figure of merit"
    % (len(files), num_failed, len(without_figure))
  )
  return 1 if num_failed else 0


if __name__ == "__main__":
  sys.exit(main())
