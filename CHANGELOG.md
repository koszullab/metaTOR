# Change Log

All notable changes to this project will be documented in this file.

## [1.3.7] - 2025-02-22
- Auto-release upon tag creation.
- Generate a Docker container hosted on `ghcr.io/koszullab/metator` upon release.

## [1.3.6] - 2025-02-21
- Prepare for bioconda release.
- Only rely on python `pypairix` package (installable from pip). 
- Binaries are now put into `./bin` instead of `./external/artifacts/`. 
- micomplete is now available from `bioconda`. 

## [1.3.4] - 2025-02-19
- Package now relies on `pyproject.toml` for build configuration with `hatch`. 
- Binaries for `louvain` and `leiden` clustering algorithms are now embedded in the package.
- Uses pre-commit hooks for code formatting and linting.
- Fix deprecated Bio.SeqUtils.GC to Bio.SeqUtils.gc_fraction.
- `Biopython` is pinned <= 1.80 to work with `micomplete 1.1.1`.

## [1.3.3] - 2023-11-27
- Improve ci.
- Add pairix as requirements.
- Correct pipeline bug.

## [1.3.2] - 2023-11-27
- Add log files for all modules.

## [1.3.1] - 2023-11-20
- Switch pairs as an argument for `metator mge`.
- Correct bug in QC.
- Add checkV in the requirements.

## [1.3.0] - 2023-11-09
- Merge MetaVIR (https://github.com/ABignaud/MetaVir) module to MetaTOR.
- Automatic generation of the compressed and sorted pairs files in the main 
pipeline.
- Correct QC issues with low-quality libraries.
- Add a conda environment file.

## [1.2.8] - 2023-04-12
- Correct issue with anvio binning file.

## [1.2.7] - 2023-04-07
- Correct bug during final plot generation.

## [1.2.6] - 2023-04-05
- Correct bug of contigs data final bin value from multiple recursive steps.

## [1.2.5] - 2023-04-05
- Remove recursive parallel which seems to failed in some case.
- Correct bug of contigs data final bin value from multiple recursive steps.
- Add exception to avoid error while mutliple metator are started from the same 
working directory.

## [1.2.4] - 2023-04-04
- Add a prefix option for fasta files.

## [1.2.3] - 2023-04-03
- Modify pairs option as an argument for QC module.
- Correct bug if only one MAG corrected in the recursive step.
- Correct figures to work with miComplete output.
- Improve quality check metrics based on mock community results.
- Update QC tests.

## [1.2.2] - 2023-03-29
- Repair the pipeline module.

## [1.2.1] - 2023-03-29
- Make scaffolding optional in pipeline.
- Remove the option to start with network file.
- Allow parallel scaffolding at the end of the pipeline.
- Force pairix index used in the contact map module.
- Change pairs from option to arguments to handle easier multiple files as 
arguments.
- Debug scaffold module default parameters.
- Update docs.

## [1.2.0] - 2023-03-17
- Remove the skip validation option.
- Use miComplete instead of checkM for the recursive step.
- Make multiple recursive step if necessary. 
- Sort pairs and do the scaffolding in the pipeline module.
- Add QC, pairs and scaffold modules.
- Correct pairs indexed loading.
- Add option to remove reads mapping at the edges of the contigs.
- Add possibility to used compressed pairs as input.
- Parallelization of the recursive step.
- Change bin id to make them sortable.
- Correct reads extraction bugs for contact map module.
- Correct bugs from cutsite mapping mode.
- Remove Looseversion sorting Warning.
- Update requirements file.
- Change output plot colors.
- Add tests for QC module.
- Update tests of network, partition and validation modules.
- Update documentation.

## [1.1.6] - 2022-10-06
- Add iterative and cutsite alignment mode.

## [1.1.5] - 2022-07-05
- Add possibility to start with bam files.
- Add SCN normalization to balance the heatmap.
- Add conda environnement and docs.
- Add prefix name for metator contact map.

## [1.1.4] - 2021-07-23

- Correct broken release 1.1.3

## [1.1.3] - 2021-07-23

- Avoid raising error if no pairs extracted in contact map module.
- Add network heatmap plot in the ouput.

## [1.1.2] - 2021-07-13

- Add the possibility to read pairix index files in the contact map module to retrieve faster pairs from a 2D region.

## [1.1.1] - 2021-06-30

Correct some minor issues:

- Correct pairs files:
  - Correct position index (previously 0-based but should be 1-based).
  - Change pairs entry to have an upper triangle shape.
- Correct recursive bin issue in the contig data file were unbinned contigs or binned contigs in bins smaller than the size threshold given were given an id of an existing bin.
- Correct issue in the contact map were contigs needed to be size sorted if a size threshold were given.
- Add some multiple network files in the output if multiple input samples were given.

## [1.1.0] - 2021-05-26

- Add some Figures in the output:
  - Completion/contamination distribution of the MAGs.
  - Distribution of the quality of the MAGs in the assembly.
  - Distribution of the GC and coverage inside each MAG.
- Add test for the new figure module.
- Add a module contactmap to build HiC matrix from some metaTOR objects such as contigs, bins, or other objects with contigs labels from the user. It could be use with metator contactmap [arguments]. The output could be use then with hicstuff view or cooler show to visualize matrix or with instagraal to try to scaffold MAGs.
- Modify the pairs files:
  - Change pairs file extension to .pairs.
  - Reorder the columns of .pairs file: readID-chr1-pos1-chr2-pos2-strand1-strand2 instead of readID-chr1-pos1-strand1-chr2-pos2-strand2.
  - Use tabulation to separate column instead of space to match the official specification.
- Correction of some minor issues:
  - Stop the duplication of the logging message.
  -Silence the HTSlib index error (from pysam.AlignementFile) which was meaningless.
  - Remove pyfastx index.
  - Put as optional the output of the clustering matrix which needs high memory usage for large assembly and optimize memory usage to build it.
  - Modify Assignment method to dataframe to avoid assignment to copy of the dataframe (SettingWithCopyWarning).

## [1.0.4] - 2021-05-14

- Change the format of the pairs file to match with the official specification.
- Add the HiC coverage (mean number of contacts in 1kb region) in the final bin summary.
- Add possibility to start with the pairs file.
- Add info about general alignment if multiples files are given.
- Organize temporary folder.

## [1.0.3] - 2021-05-04
  
- Add clustering matrix in the output of the partition, validation and pipeline modules.
- Add output files documentation.
- Add codeco.

## [1.0.2] - 2021-05-01

- Modification of the ouput files:
  - Add bin coverage based on the depth file given in the bin_summary file.
  - Add a final bin column in the contig data file after the validation to specify the name of the final bin of the contig. A "ND" value is given if the contig isn't binned.
  - Add the pairs file of the alignment to be able to produce HiC matrices from it using hicstuff for example.
  - Add a log file with all information of the run.
  - With the pipeline workflow remove intermediate contig data files.
- Modify call to input data to use only header names and not order of the column.
- Change the parameter of the bowtie2 alignement back to --very-sensitive-local.
- Add the total size of the different category of bins in the log.
- Correct minor issues in the doc.

## [1.0.1] - 2021-04-27

Correct bug in validation module:

- louvain_recursif: correct type error to identify contigs of a subnetwork which could triggers empty subnetwork error.

## [1.0.0] - 2021-04-26

New version of metaTOR in Python3
