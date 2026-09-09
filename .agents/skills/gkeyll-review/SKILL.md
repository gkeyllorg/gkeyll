---
name: gkeyll-review
description: Principles of writing good code, which would pass review. Agents must review their own code.
---

# Gkeyll code review principles

Think carefully and critically, each point one at a time, whether the code written adheres to the standards of the project.

- There are no if statements inside time loops. Instead, use a function pointer to select the appropriate function.
- Code must be formatted with clang-format.
