---
name: naming-conventions
description: Apply Gkeyll naming conventions when creating, editing, or reviewing C, CUDA, or Lua files and code elements.
---

# Instructions

* Apply these naming conventions when creating, editing, or reviewing C, CUDA, and Lua code.

# Naming conventions

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

- Do not use single letter names for variables whose scope spans more than 15 lines.
