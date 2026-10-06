---
name: software-design
description: Some design practices in Gkeyll code to consider when implementing new modules (apps, updaters) or modifying existing ones.
---

# Instructions

* Load this skill whenever new modules or updaters in zero/ or apps/ folders
  are to be created, or when existing ones will be edited.

## Software design elements

### Module structure

Most modules in zero/ or apps/ consist of 3 public functions:
1. A creation function, typically called gkyl_<module>_new or gkyl_<module>_init.
2. An execution function, often called gkyl_<module>_advance or gkyl_<module>_apply.
3. A deletion function, typically called gkyl_<module>_release.

There may also be some additional auxiliary private or public functions.

### Module best practices

- All dynamic (heap) allocations should happen in the initialization function,
  and freed in the release function, using the appropriate gkyl_ allocation/deallocation functions when possible.
- Don't place logic branching (e.g. if-statements) that depend on a
  time-independent choice or parameter inside the methods called in the time loop of a simulation (e.g. _advance). Instead, use function pointers to set the appropriate method during the initialization of the module, and call that method inside the time loop.

### Other principles to follow

- Consider extensibility, maintainability, simplicity and how modular design.
- Avoid code duplication whenever possible (e.g. write functions called
  multiple times) and without breaking layering.
- Write shorter code and refactor into a sub-module whenever possible.
