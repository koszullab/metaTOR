"""Run miComplete with patches for recent Biopython and Python, and with
pyrodigal and pyhmmer in place of the prodigal and hmmsearch programs.

miComplete (<= 1.1.1, no longer updated upstream) relies on two names that
have since been removed:

    - ``Bio.SeqUtils.GC``, removed in Biopython 1.82 (replaced by
      ``Bio.SeqUtils.gc_fraction``);
    - ``distutils.spawn``, removed from the standard library in Python 3.12.

This module puts both names back before miComplete is imported, so that
metator does not have to pin an old Biopython.

miComplete also runs two programs that pip does not install: prodigal, which
finds the genes of each bin, and hmmsearch (HMMER), which finds the marker
genes among their proteins. This module has miComplete use pyrodigal and
pyhmmer instead, the Python bindings to Prodigal and HMMER, which pip installs
with metator: the bins are validated without prodigal or HMMER installed.

The patches are applied at import time (not only under ``__main__``) so that
the worker processes of miComplete's multiprocessing pool also get them
whatever the start method. They are only applied in the miComplete
subprocess: the Python process running metator is left untouched.

Usage::

    python -m metator._micomplete <miComplete arguments>
"""

import functools
import logging
import os
import shutil
import sys
import types
import warnings

import pyhmmer
import pyrodigal
from Bio import SeqIO

# The programs that miComplete runs, and the modules that stand in for them.
PROGRAMS = {"prodigal": pyrodigal, "hmmsearch": pyhmmer}


def patch_biopython_gc():
    """Restore ``Bio.SeqUtils.GC`` if the installed Biopython lacks it."""
    import Bio.SeqUtils

    if hasattr(Bio.SeqUtils, "GC"):
        return

    def GC(seq):
        """GC content in percent, computed as ``Bio.SeqUtils.GC`` did in
        Biopython <= 1.81 (this is its body): G, C and S are counted in either
        case, every character (including ambiguous ones) counts in the length,
        and an empty sequence gives 0.0.

        ``100 * gc_fraction(seq, ambiguous="ignore")`` counts the same
        characters but divides before multiplying, which can differ in the
        last bit and so move miComplete's GC-content, rounded to two decimals,
        by 0.01.
        """
        gc = sum(seq.count(x) for x in ["G", "C", "g", "c", "S", "s"])
        try:
            return gc * 100.0 / len(seq)
        except ZeroDivisionError:
            return 0.0

    Bio.SeqUtils.GC = GC


def patch_distutils_spawn():
    """Provide ``distutils.spawn.find_executable`` if distutils is missing."""
    try:
        from distutils import spawn
    except ImportError:

        def find_executable(executable, path=None):
            return shutil.which(executable, path=path)

        spawn = types.ModuleType("distutils.spawn")
        spawn.find_executable = find_executable
        distutils = types.ModuleType("distutils")
        distutils.spawn = spawn
        sys.modules["distutils"] = distutils
        sys.modules["distutils.spawn"] = spawn


def create_proteome(fasta, base_name=None, logger=None):
    """Find the genes of a bin with pyrodigal, where miComplete runs
    ``prodigal -i <fasta> -a <base_name>_prodigal.faa``.

    As prodigal in its single genome mode, its default: trained on all the
    contigs of the fasta file together, then run on each contig. The proteins
    of the genes are written as prodigal writes them, to
    <base_name>_prodigal.faa in the current directory. With less than 20 kb of
    sequence, too little to train on, the file is left empty, as prodigal
    leaves it: miComplete then finds no gene and no marker in the bin.

    Parameters
    ----------
    fasta : str
        Path to the fasta file of a bin.
    base_name : str
        Name of the proteome file, before "_prodigal.faa". [Default: the name
        of the fasta file without its extension]
    logger : logging.Logger
        Unused, as in miComplete.

    Returns
    -------
    str:
        Path to the fasta file of the proteins.
    """
    if not base_name:
        base_name = "".join(os.path.basename(fasta).split(".")[:-1])
    prot_filename = base_name + "_prodigal.faa"
    contigs = [(record.id, bytes(record.seq)) for record in SeqIO.parse(fasta, "fasta")]
    with open(prot_filename, "w") as proteome:
        if not contigs:
            return prot_filename
        finder = pyrodigal.GeneFinder(meta=False)
        try:
            with warnings.catch_warnings():
                # Under 100 kb, as bins often are. prodigal warns too, on a
                # stderr that miComplete discards.
                warnings.filterwarnings("ignore", message="sequence should be at least", category=UserWarning)
                # The contigs are joined as prodigal joins them.
                finder.train(*(seq for _, seq in contigs))
        except ValueError:
            # Under 20 kb.
            return prot_filename
        for name, seq in contigs:
            finder.find_genes(seq).write_translations(proteome, sequence_id=name)
    return prot_filename


@functools.cache
def read_hmms(path):
    """The profile HMMs of an HMMER file, read once per process."""
    with pyhmmer.plan7.HMMFile(path) as hmm_file:
        return tuple(hmm_file)


def hmm_search(self):
    """Search the marker HMMs of miComplete in the proteins of a bin with
    pyhmmer, where ``micomplete.completeness.calcCompleteness.hmm_search``
    runs ``hmmsearch -E <evalue> --tblout <tblout> <hmms> <proteome>``.

    The hits are written as hmmsearch writes them, in the --tblout table that
    miComplete reads. A proteome without any protein (prodigal found no gene)
    is refused, as hmmsearch refuses it: with a warning, no table, and 1
    returned rather than 0.

    Returns
    -------
    str:
        Path to the table of the hits.
    int:
        0, or 1 when the proteome holds no protein.
    """
    try:
        self.logger.log(logging.INFO, "Starting hmmsearch")
    except AttributeError:
        pass
    hmms = read_hmms(self.hmms)
    try:
        with pyhmmer.easel.SequenceFile(self.fasta, format="fasta", digital=True, alphabet=hmms[0].alphabet) as proteome:
            proteins = proteome.read_block()
    except EOFError:
        message = f"Error thrown by HMMER, is {self.fasta} empty?"
        try:
            self.logger.log(logging.WARNING, message)
        except AttributeError:
            print(f"Warning: {message}", file=sys.stderr)
        return self.tblout, 1
    # calcCompleteness keeps the e-value as the hmmsearch option: "-E <evalue>".
    evalue = float(self.evalue.split()[-1])
    with open(self.tblout, "wb") as tblout:
        # One CPU: miComplete already runs one process per thread.
        for i, hits in enumerate(pyhmmer.hmmsearch(hmms, proteins, cpus=1, E=evalue)):
            hits.write(tblout, format="targets", header=i == 0)
    return self.tblout, 0


def which(cmd, mode=os.F_OK | os.X_OK, path=None):
    """``shutil.which`` for miComplete, which checks that prodigal and
    hmmsearch are on the PATH before it starts: they are found, installed or
    not, as the modules that stand in for them."""
    if cmd in PROGRAMS:
        return PROGRAMS[cmd].__file__
    return shutil.which(cmd, mode=mode, path=path)


def patch_prodigal_and_hmmsearch():
    """Have miComplete run pyrodigal and pyhmmer in place of the prodigal and
    hmmsearch programs. This imports miComplete: patch_biopython_gc and
    patch_distutils_spawn must come first."""
    import micomplete.micomplete
    from micomplete.completeness import calcCompleteness

    cli = micomplete.micomplete
    cli.create_proteome = create_proteome  # ty: ignore[invalid-assignment]
    calcCompleteness.hmm_search = hmm_search
    # A copy of shutil, with the which() above, for miComplete only.
    cli_shutil = types.ModuleType("shutil")
    cli_shutil.__dict__.update(vars(shutil))
    cli_shutil.which = which  # ty: ignore[missing-attribute]
    cli.shutil = cli_shutil  # ty: ignore[invalid-assignment]


patch_biopython_gc()
patch_distutils_spawn()
patch_prodigal_and_hmmsearch()


def main():
    """Run the miComplete command line interface."""
    from micomplete.micomplete import main as micomplete_main

    sys.argv[0] = "miComplete"
    return micomplete_main()


if __name__ == "__main__":
    sys.exit(main())
