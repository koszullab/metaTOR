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

## Versioning

- We follow [conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) for versioning this project.
- We use [semantic versioning](https://semver.org) for this project.
- Whenever a version bump is needed, the maintainer will bump the version in [pyproject.toml](./pyproject.toml), and tag the commit. This will trigger: 
  - A new release on Github; 
  - Propagation of this release on PyPI; 
  - A new Docker image on [ghcr.io/koszullab/metator](https://ghcr.io/koszullab/metator).

## Acknowledgments

This document is based off of the [guidelines from the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).

