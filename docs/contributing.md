# Contributing

Bug reports, suggestions, documentation and pull requests are welcome on
[GitHub](https://github.com/koszullab/metaTOR). We follow the
[fork and pull model](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models),
[conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) and
[semantic versioning](https://semver.org).

## Development set-up

The environment is managed by [uv](https://docs.astral.sh/uv/); `uv.lock`
pins every dependency. Install uv, then, from a clone of your fork:

```bash
git clone https://github.com/koszullab/metaTOR.git
cd metaTOR
uv sync                       # Python 3.13 (.python-version), the dev tools, metator in editable mode
uv run pre-commit install     # optional: ruff on every commit
```

`uv run metator --version` runs the version in the repository.

## Tests, lint and format

```bash
uv run pytest                     # the tests, with coverage
uv run ruff check src tests       # lint
uv run ruff format src tests      # format (line length 130)
uv run ty check src               # type check (advisory for now)
```

- New functions come with tests in `tests/`. Codecov requires at least 80%
  coverage of the lines a pull request changes.
- [CI](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/ci.yml)
  runs on every pull request and push to `master`: the lock file, ruff, ty,
  the tests on each supported Python version, and a job on the newest
  releases of the dependencies (also weekly).
- Docstrings follow the [numpydoc standard](https://numpydoc.readthedocs.io/en/latest/format.html):
  they make the [Python API reference](api/network.md).

## Documentation

The site is built from `docs/` by [Zensical](https://zensical.org) and served
at <https://koszullab.github.io/metator/>.

```bash
uv sync --extra docs
uv run zensical serve     # live preview at http://localhost:8000
uv run zensical build     # the site, in site/
```

- **Pages**: the Markdown files of `docs/`, listed in order by `nav` in
  `zensical.toml`. A new page goes in both.
- **Python API reference**: the pages of `docs/api/` hold
  `::: metator.<module>` directives, which
  [mkdocstrings](https://mkdocstrings.github.io) expands from the docstrings.
- **CLI reference**: `docs/cli.md`, written by hand. `tests/test_docs.py`
  fails when it misses a command or option of the CLI, or documents one the
  CLI lacks.
- **Diagrams**: [Mermaid](https://mermaid.js.org) code blocks. **Images**:
  `docs/images/`.
- **Look**: [chromokit](https://github.com/js2264/chromokit), the UI kit of
  metacooler's documentation. Do not edit its files (`overrides/`,
  `docs/stylesheets/chromokit.css`, `docs/javascripts/chromokit.js`, listed
  in `chromokit.json`); put metaTOR's own styles in
  `docs/stylesheets/site.css`. To update the kit, and check that no kit file
  was edited:

    ```bash
    npx github:js2264/chromokit sync
    npx github:js2264/chromokit check
    ```

Every push to `master` builds and deploys the site
([docs.yml](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/docs.yml));
pull requests build it without deploying. The build is strict: any warning,
such as a broken link, fails it.

## Changelog

Every user-visible change gets an entry under `## Unreleased` in
`docs/changelog.md`, in *Added*, *Changed*, *Removed* or *Fixed*.

## Releasing

Releases use [bump-my-version](https://github.com/callowayproject/bump-my-version),
configured in `pyproject.toml`. From an up-to-date `master` with green CI:

```bash
uv run bump-my-version bump patch   # or minor / major
git push --follow-tags
```

The bump sets `__version__` in `src/metator/version.py` (the single source of
the version) and `release` in `zensical.toml` (the version the site shows),
turns `## Unreleased` in `docs/changelog.md` into a dated section, commits as
`chore(release): vX.Y.Z` and tags `vX.Y.Z`. Preview it with
`uv run bump-my-version bump patch --dry-run -vv`.

Pushing the tag starts
[release.yml](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/release.yml),
which:

- checks that the tag matches `__version__`, builds the sdist and wheel, and
  checks that `metator` runs from the wheel in a clean environment;
- publishes them to PyPI through
  [trusted publishing](https://docs.pypi.org/trusted-publishers/);
- creates the GitHub release from the version's section of
  `docs/changelog.md`;
- builds the Docker image from the wheel, checks it, and pushes it to
  [ghcr.io/koszullab/metator](https://ghcr.io/koszullab/metator).

The release commit on `master` deploys the site with the new version.
Bioconda picks up the PyPI release itself, but its recipe lists the run
dependencies by hand: update them in
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/metator)
when they change.

This guide is based on the
[guidelines of the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).
