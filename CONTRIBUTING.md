# Contributing to de_quack

Thanks for your interest in contributing! This document covers how to
set up a development environment, run tests, and submit changes.

## Current maintenance status

de_quack is currently maintained by a single developer. Response times
to issues and PRs may vary. Contributions, bug reports, and feature
requests are all welcome.

## Setting up a development environment

1. Clone the repository:
```bash
   git clone https://github.com/yangp7833-commits/De_quack.git
   cd De_quack
```
2. Install in editable mode with development dependencies:
```bash
   pip install pytest
   pip install -e
```
(A `[dev]` extra will be added here once one exists.)

## Running tests

```bash
pytest
```

## Reporting bugs

Open a GitHub issue with:
- A minimal example that reproduces the problem
- The de_quack version (`pip show de_quack`) and Python version
- The full error traceback, if applicable

## Suggesting features

Open an issue describing the use case, not just the feature — it
helps to know what problem you're trying to solve.

## Submitting changes

1. Fork the repo and create a branch from `main`.
2. Make your changes, with tests for any new behavior.
3. Make sure `pytest` passes locally.
4. Open a pull request with a clear description of what changed and why.

## Code style

No code style is enforced at this time.

## Questions

Feel free to open an issue for questions too — no need for it to be
a bug or a concrete feature request.
