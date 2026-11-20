# CLAUDE.md

Guidance for Claude Code (claude.ai/code) working in this repository.

## Pull requests

The maintainer squash-merges every pull request, so a pull request is one
change, whatever its commits:

- **Reviews:** when asked to review a pull request, review it as a whole, not
  commit by commit: say how the whole change affects the package (its
  behaviour, outputs, command line, Python API and docs). This applies to the
  brief given to an independent reviewer too.
- **Tests:** run the test suite on the head of the pull request only, not at
  every commit.
