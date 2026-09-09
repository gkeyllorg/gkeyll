---
name: gkeyll-naming-conventions
description: Apply Gkeyll file, function, and variable naming conventions and module lifecycle patterns when creating, editing, or reviewing C, CUDA, and Lua code.
---

# Gkeyll naming conventions

### Files.

- Public header files have names starting with gkyl_ and ending in .h.
- Private header files have names starting with gkyl_ and ending in _priv.h.
- CUDA files have names ending in _cu.cu. 
- luareg/ folders have Lua input files whose names end with .lua. The Lua
  wrappers are in the apps/ folders and have names ending in _lw.c

### Functions

- Public functions (defined in public header files) should have a name that starts
  with gkyl_.
- Public functions in files in zero/ folders should have a name that starts
  with the name of the file that contains it.
- Private functions (defined in private header files) or static (and not defined
  in private headers) in files in zero/ folders should have a name that starts
  with an abbreviated version of the name of the file that contains it.

### Variables

- Do not use single letter names for variables whose scope spans more than ~15
  lines.

## Module structure and pattern

Most modules in zero/ or apps/ consist of 3 public functions:
1. A creation function, typically called gkyl_<module>_new or gkyl_<module>_init.
2. An execution function, typically called gkyl_<module>_advance or gkyl_<module>_apply.
3. A deletion function, typically called gkyl_<module>_release.

There may also be some additional auxiliary private or public functions.

## Modular design

Each file must have a clear descriptive name.
Filenames help organize the code into the relevant sections.
When a set of code could be made smaller and refactored into a sub-module, it should.
The purpose and intent of the file should be clear from its name.