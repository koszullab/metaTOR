#!/usr/bin/env python3

"""Former entry point of the metator command, kept as an alias of
:func:`metator.cli.cli`: ``python -m metator.main`` and ``metator.main.main()``
still run the command line.
"""

from metator.cli import cli

main = cli

if __name__ == "__main__":
    main()
