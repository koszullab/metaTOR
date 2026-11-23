# Contributing

All contributions are welcome: bug reports, suggestions, documentation and
pull requests, on [GitHub](https://github.com/koszullab/metaTOR). We follow
the [fork and pull model](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models),
[conventional commits](https://www.conventionalcommits.org/en/v1.0.0/) and
[semantic versioning](https://semver.org).

## Development set-up

metaTOR's development environment is managed by
[uv](https://docs.astral.sh/uv/), and `uv.lock` pins every dependency. Install
uv; then, from a clone of your fork:

```bash
git clone https://github.com/koszullab/metaTOR.git
cd metaTOR
uv sync                       # Python 3.13 (.python-version), the dev tools, metator in editable mode
uv run pre-commit install     # optional: ruff on every commit
```

`uv run metator --version` then runs the version in the repository.

## Tests, lint and format

```bash
uv run pytest                     # the tests, with coverage
uv run ruff check src tests       # lint
uv run ruff format src tests      # format (line length 130)
uv run ty check src               # type check (advisory for now)
```

New functions should come with tests in `tests/`, and the lines a pull request
changes must be at least 80% covered: Codecov checks it.
[CI](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/ci.yml)
runs on every pull request and every push to `master`: the lock file, ruff,
ty, the tests on each Python version metaTOR supports, and a job that
re-resolves the dependencies to their newest releases (also weekly).

Docstrings follow the [numpydoc standard](https://numpydoc.readthedocs.io/en/latest/format.html):
they are the [Python API reference](api/network.md) of this site.

## Documentation

This site is built from the Markdown files of `docs/` by
[Zensical](https://zensical.org), and served by GitHub Pages at
<https://koszullab.github.io/metator/>.

```bash
uv sync --extra docs
uv run zensical serve     # live preview at http://localhost:8000
uv run zensical build     # the site, in site/
```

- **Pages** are the Markdown files of `docs/`; `nav` in `zensical.toml` lists
  them in order. A new page goes in both.
- **The Python API reference** is written from the docstrings: each page of
  `docs/api/` holds `::: metator.<module>` directives, which
  [mkdocstrings](https://mkdocstrings.github.io) expands.
- **The CLI reference**, `docs/cli.md`, is written by hand. `tests/test_docs.py`
  fails when a command or an option of the CLI is missing from it, or when it
  documents one the CLI does not have.
- **Diagrams** are [Mermaid](https://mermaid.js.org) code blocks.
- **Images** are in `docs/images/`.
- **The look** comes from [chromokit](https://github.com/js2264/chromokit),
  the UI kit of metacooler's documentation. Its files (`overrides/`,
  `docs/stylesheets/chromokit.css`, `docs/javascripts/chromokit.js`) are
  listed in `chromokit.json`: do not edit them, put metaTOR's own styles in
  `docs/stylesheets/site.css`. To take up a newer chromokit, and to check that
  no kit file was edited:

    ```bash
    npx github:js2264/chromokit sync
    npx github:js2264/chromokit check
    ```

Every push to `master` builds the site and deploys it
([docs.yml](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/docs.yml));
pull requests to `master` build it, without deploying. The build is strict:
a warning, such as a link to a page that does not exist, fails it.

## Changelog

Every change that a user could notice gets an entry under `## Unreleased` in
`docs/changelog.md`, in its section: *Added*, *Changed*, *Removed* or
*Fixed*.

## Releasing

Releases are cut with [bump-my-version](https://github.com/callowayproject/bump-my-version),
configured in `pyproject.toml`. From an up-to-date `master` whose CI is green:

```bash
uv run bump-my-version bump patch   # or minor / major
git push --follow-tags
```

The bump sets `__version__` in `src/metator/version.py`, the single source of
the version, and `release` in `zensical.toml`, the version this site shows;
turns `## Unreleased` in `docs/changelog.md` into a dated section for the
version; commits as `chore(release): vX.Y.Z` and creates the `vX.Y.Z` tag.
Preview it with `uv run bump-my-version bump patch --dry-run -vv`.

Pushing the tag starts
[release.yml](https://github.com/koszullab/metaTOR/blob/master/.github/workflows/release.yml),
which:

- checks that the tag matches `__version__`, builds the sdist and the wheel,
  and installs the wheel in a clean environment to check that `metator` runs;
- publishes them to PyPI through
  [trusted publishing](https://docs.pypi.org/trusted-publishers/), with no
  stored token;
- creates the GitHub release, with the section of `docs/changelog.md` for the
  version;
- builds the Docker image from the released wheel, checks that it runs, and
  pushes it to [ghcr.io/koszullab/metator](https://ghcr.io/koszullab/metator).

The push of the release commit to `master` deploys this site, with the new
version in its header. Bioconda picks up the new PyPI release by itself, but
its recipe lists the run dependencies by hand: when they change, update them
in [bioconda-recipes](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/metator).

This guide is based on the
[guidelines of the sparse project](https://github.com/pydata/sparse/blob/master/docs/contributing.rst).
