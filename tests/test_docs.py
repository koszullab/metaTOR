"""The CLI reference of the documentation (docs/cli.md) follows the CLI."""

import re
from pathlib import Path

import click
import pytest

from metator.cli import cli

CLI_MD = Path(__file__).resolve().parents[1] / "docs" / "cli.md"

#: The heading of a command's section: ### `pipeline` { #pipeline }
_SECTION = re.compile(r"^### `([a-z0-9_-]+)`", re.MULTILINE)
#: The first cell of an options table row: | `-a, --assembly` | ...
_OPTION_CELL = re.compile(r"^\| `(-[^`]+)` \|", re.MULTILINE)


def documented_options():
    """The options documented for each command, by command name."""
    text = CLI_MD.read_text(encoding="utf-8")
    headings = list(_SECTION.finditer(text))
    sections = {}
    for heading, following in zip(headings, [*headings[1:], None], strict=True):
        body = text[heading.end() : following.start() if following else len(text)]
        options = set()
        for cell in _OPTION_CELL.findall(body):
            options.update(name.strip() for name in cell.split(","))
        sections[heading.group(1)] = options
    return sections


def cli_options(command):
    """The options of a click command, without --help."""
    options = set()
    for param in command.params:
        if isinstance(param, click.Option):
            options.update(param.opts + param.secondary_opts)
    return options - {"-h", "--help"}


def test_every_command_is_documented():
    ctx = click.Context(cli)
    assert set(documented_options()) == set(cli.list_commands(ctx))


@pytest.mark.parametrize("name", sorted(cli.list_commands(click.Context(cli))))
def test_every_option_is_documented(name):
    command = cli.get_command(click.Context(cli), name)
    documented = documented_options()[name]
    missing = cli_options(command) - documented
    unknown = documented - cli_options(command)
    assert not missing, f"options of `metator {name}` missing from docs/cli.md: {missing}"
    assert not unknown, f"docs/cli.md documents options `metator {name}` lacks: {unknown}"
