# Contributing

## General guidelines

We follow the [fork and pull model](https://help.github.com/articles/about-collaborative-development-models) on GitHub. We also recommend checking out this recommended [git workflow](https://www.asmeurer.com/git-workflow/).

## Contributing Code

This project has a number of requirements for all code contributed.

* We use [uv](https://docs.astral.sh/uv/) to manage the development environment; `uv.lock` pins every dependency.
* We use [pytest](https://docs.pytest.org/en/latest) for unit tests.
* We use [ruff](https://docs.astral.sh/ruff/) for linting and code formatting.
* We use [ty](https://docs.astral.sh/ty/) for static type checking (advisory for now).

[CI](.github/workflows/ci.yml) runs on every push and pull request: lock-file check, ruff, ty, and the test suite on each supported Python version from `uv.lock`, plus a job that re-resolves to the newest dependency releases (also run weekly).

## Setting up Your Development Environment

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), and the non-Python tools metator calls (uv cannot install them):

```sh
sudo apt install bowtie2 bwa samtools hmmer prodigal   # or: brew install bowtie2 bwa samtools hmmer prodigal
```

Then, after forking and cloning the repository:

```sh
git clone https://github.com/koszullab/metaTOR.git
cd metaTOR
uv sync                       # Python 3.11 (.python-version) + dev tools, metator in editable mode
uv run pre-commit install     # optional: ruff on every commit
```

## Running/Adding Unit Tests

From the root of the repository:

```sh
uv run pytest                     # tests, with coverage
uv run ruff check src tests       # lint
uv run ruff format src tests      # format
uv run ty check src               # type check
```

Unit tests are automatically run on GitHub for all pull requests.

## Pull Requests

Pull requests are welcome, once changes have been tested and pass all CI checks. Please make sure to include tests for any new functionality. 

## Releasing

- We follow [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) and [semantic versioning](https://semver.org).
- Changelog entries for unreleased changes go under `## Unreleased` in [CHANGELOG.md](./CHANGELOG.md) as they are merged.
- To cut a release, from an up-to-date `master` whose CI is green:

  ```sh
  uv run bump-my-version bump patch   # or minor / major
  git push --follow-tags
  ```

  [bump-my-version](https://github.com/callowayproject/bump-my-version) (configured in [pyproject.toml](./pyproject.toml)) promotes `## Unreleased` to a dated section, bumps `__version__` in `src/metator/version.py` (the single source of truth for the version), then commits and tags. Pushing the tag triggers [release.yml](.github/workflows/release.yml), which:
  - checks that the tag matches `__version__`, builds the sdist and wheel, and installs the wheel in a clean environment to check that `metator` runs;
  - publishes them to PyPI through [trusted publishing](https://docs.pypi.org/trusted-publishers/) (no stored token);
  - creates the GitHub release with this version's changelog section;
  - builds the Docker image from the released wheel, checks that it runs and pushes it to [ghcr.io/koszullab/metator](https://ghcr.io/koszullab/metator).
- Bioconda picks up the new PyPI release automatically, but its recipe lists the run dependencies by hand: when they change, update them in [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/metator).

## Acknowledgments

This document is based off of the [guidelines from the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).

