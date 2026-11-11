# Contributing

## General guidelines

We follow the [fork and pull model](https://help.github.com/articles/about-collaborative-development-models) on GitHub. We also recommend checking out this recommended [git workflow](https://www.asmeurer.com/git-workflow/).

## Contributing Code

This project has a number of requirements for all code contributed.

* We use [pytest](https://docs.pytest.org/en/latest) for unit tests.
* We use [Black](https://black.readthedocs.io/en/stable/) for code formatting.
* We use [mypy](http://mypy-lang.org) for static type checking. 

We strongly enforce continuous integration (CI) on all pull requests. This includes running unit tests, code coverage, code style checks, and static type checks. We use [Github Actions](.github/workflows/ci.yml) for this purpose. 

## Setting up Your Development Environment

After forking and cloning the repository, install in "editable" (i.e. development) mode using the `-e` option:

```sh
git clone https://github.com/koszullab/metator.git
cd metator
pip install -e .[dev]
```

## Running/Adding Unit Tests

We use [pytest](https://docs.pytest.org/en/latest) as our unit testing framework. Once you've configured your environment, you can just `cd` to the root of your repository and run

```sh
pytest
```

Unit tests are automatically run on GitHub for all pull requests.

## Pull Requests

Pull requests are welcome, once changes have been tested and pass all CI checks. Please make sure to include tests for any new functionality. 

## Releasing

- We follow [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) and [semantic versioning](https://semver.org).
- Changelog entries for unreleased changes go under `## Unreleased` in [CHANGELOG.md](./CHANGELOG.md) as they are merged.
- To cut a release, from an up-to-date `master` whose CI is green:

  ```sh
  uvx bump-my-version bump patch   # or minor / major
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

