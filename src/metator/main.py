#! /usr/bin/env python3
# Based on Rémy Greinhofer (rgreinho) tutorial on subcommands in docopt
# https://github.com/rgreinho/docopt-subcommands-example
# abignaud, 20201118

"""
MetaHiC pipeline for generating and manipulating MAGs.

usage:
    metator [-hv] <command> [<args>...]

options:
    -h, --help              shows the help
    -v, --version           shows the version

The subcommands are:
    network     Generate metaHiC contigs network from fastq reads or bam files
                and normalize it.
    partition   Partition metaHiC network using Louvain or Leiden algorithm.
    pipeline    Use all the others command to give a binned output in one
                command line.
    validation  Validates bins using CheckM and make a recursive partition to
                try to decontaminate them.
    qc          Generates some quality check on the output of metator.
    contactmap  Generates a HiC contact map from one metaTOR object from the
                final ouptut of metaTOR.
    scaffold    Scaffold a metator bin based on pairs files.
    pairs       Sort the pairs file using pairtools. Compress them using bgzip.
                Index them using pypairix.
    host        Detect bacterial host from a metaHiC network binned by metaTOR
                given a annotated MGE list.
    mge         Build MGE MAGs based on metagenomic binning using metabat2
                or pairs from the metaHiC data.
"""

from docopt import docopt
from docopt import DocoptExit
import metator.commands as commands
from metator.version import __version__


def main():
    args = docopt(__doc__, version=__version__, options_first=True)
    # Retrieve the command to execute.
    command_name = args.pop("<command>").capitalize()

    # Retrieve the command arguments.
    command_args = args.pop("<args>")
    if command_args is None:
        command_args = {}
    # After 'popping' '<command>' and '<args>', what is left in the args
    # dictionary are the global arguments.

    # Retrieve the class from the 'commands' module.
    try:
        command_class = getattr(commands, command_name)
    except AttributeError:
        print("Unknown command.")
        raise DocoptExit()
    # Create an instance of the command.
    command = command_class(command_args, args)
    # Execute the command.
    command.execute()


if __name__ == "__main__":
    main()
