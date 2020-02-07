[![Build Status](https://travis-ci.org/fulcrumgenomics/fgbio.svg?branch=master)](https://travis-ci.org/fulcrumgenomics/fgbio)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/fgbio.svg?label=Bioconda)](http://bioconda.github.io/recipes/fgbio/README.html)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/nebiolabs/tasmanian/blob/master/LICENSE)
[![Language](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org)

<!-- Cartoon of tasmanian devil -->
<!--[](https://upload.wikimedia.org/wikipedia/en/c/c4/Taz-Looney_Tunes.svg) -->
<img src="https://github.com/nebiolabs/tasmanian/blob/master/figures/kissclipart-cartoon-devil-bugs-bunny-clipart-tasmanian-devil-b-d7d7781bbbc9a4a9.png" width="100" height="100">
<!-- image from https://www.kissclipart.com/cartoon-devil-bugs-bunny-clipart-tasmanian-devil-b-ryiinn/ NO COPYRIGHTS-->


Tasmanian
====

A tool for the analysis of artifacts (mutations with a systematic stamp/pattern) in high throughput sequencing data from DNA samples. This is specially relevant in allele variations studies.

<!---toc start-->
  * [Goals](#goals)
  * [Overview](#overview)
  * [List of tools](#list-of-tools)
  * [Building](#building)
  * [Command line](#command-line)
  * [Include fgbio in your project](#include-fgbio-in-your-project)
  * [Contributing](#contributing)
  * [Authors](#authors)
  * [License](#license)

<!---toc end-->

# Goals

The main goal is to identify systematic missmatches that might confound SNPs or other variations that should or should not be associated to biological outcomes. Since we noticed a set of regions, which might not necessarily be missplaced in the reference genome, have dramatic effects in this analysis, we provide a way of spliting these reads and incorporate the information in different tables, so that intersecting/non intersecting reads are not filtered out. Also, the researcher has a more accurate picture of the influence of these regions in the observed artifacts.

* Easy to use command-line and next-flow implementation.
* Includes a Galaxy wrapper 
* Clear and thorough documentation for each tool.
* Open source development for the benefit of the community and our clients.

## Overview

Fgbio is a set of command line tools to perform bioinformatic/genomic data analysis. 
The collection of tools within `fgbio` are used by our customers and others both for ad-hoc data analysis and within production pipelines.
These tools typically operate on read-level data (ex. FASTQ, SAM, or BAM) or variant-level data (ex. VCF or BCF).
They range from simple tools to filter reads in a BAM file, to tools to compute consensus reads from reads with the same molecular index/tag.
See the [list of tools](#list-of-tools) for more detail on the tools

## List of tools

For a full list of available tools please see the [tools section](http://fulcrumgenomics.github.io/fgbio/tools/latest) of the project website.

Below we highlight a few tools that you may find useful.

* Tools for working with Unique Molecular Indexes (UMIs, aka Molecular IDs or MIDs). 
  * Annotating/Extract Umis from read-level data: `AnnotateBamWithUmis` and `ExtractUmisFromBam`.
  * Tools to manipulate read-level data containing Umis: `CorrectUmis`, `GroupReadsByUmi`, `CallMolecularConsensusReads` and `CallDuplexConsensusReads`
* Tools to manipulate read-level data:
  * FastqManipulation: `DemuxFastqs` and `FastqToBam`
  * Filter read-level data: `FilterBam`.
  * Clipping of reads: `ClipBam`.
  * Randomize the order of read-level data: `RandomizeBam`.
  * Update read-level metadata: `SetMateInformation` and `UpdateReadGroups`.
* Quality assessment tools:
  * Detailed substitution error rate evaluation: `ErrorRateByReadPosition`
  * Sample pooling QC: `EstimatePoolingFractions`
  * Splice-aware insert size QC for RNA-seq libraries: `EstimateRnaSeqInsertSize`
  * Assessment of duplex sequencing experiments: `CollectDuplexSeqMetrics`
* Miscellaneous tools:
  * Pick molecular indices (ex. sample barcodes, or molecular indexes): `PickIlluminaIndices` and `PickLongIndices`.
  * Convert the output of HAPCUT (a tool for phasing variants): `HapCutToVcf`.
  * Find technical or synthetic sequences in read-level data: `FindTechnicalReads`.
  * Assess phased variant calls: `AssessPhasing`.

## Building 
### Cloning the Repository

[Git LFS](https://git-lfs.github.com/) is used to store large files used in testing fgbio.  In order to compile and run tests it is necessary to [install git lfs](https://git-lfs.github.com/).  To retrieve the large files either:

1. Clone the repository _after_ installing git lfs, or
2. In a previously cloned repository run `git lfs pull` once

After initial setup regular git commands (e.g. `pull`, `fetch`, `push`) will also operate on large files and no special handling is needed.

To clone the repository: `git clone https://github.com/fulcrumgenomics/fgbio.git`

### Running the build
fgbio is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.13/```.

Tests may be run with ```sbt test```.

Java SE 8 is required.


## Command line

`java -jar target/scala-2.13/fgbio-<version>.jar` to see the commands supported.  Use `java -jar target/scala-2.13/fgbio-<version>.jar <command>` to see the help message for a particular command.

## Contributing

Contributions are welcome and encouraged.
We will do our best to provide an initial response to any pull request or issue within one-week.
For urgent matters, please contact us directly.

## Authors

* [Brad Langhorst](https://github.com/bwlang) (comander/strategist)
* [Ariel Erijman](https://github.com/aerijman) (worker/maintainer)

## License

`fgbio` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE).



=======



#### TODO list:
* Adapt module report_intersect to be standalone or called by tasmanian.
* write a test suite to automatically check that everything works. --> nose test

