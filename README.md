<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Language-Python3.7-green.svg)
[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer)

Authors: Martin H&ouml;lzer, Maximilian Arlt, Eva Aßmann

# HyPro
Protein-coding annotation extension using additional homology searches against larger databases.

## Summary

The HyPro tool extends common protein-coding annotations made with [Prokka](https://github.com/tseemann/prokka) using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (available are [UniProtKB](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/), [Uniref50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/), [Uniref90](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/), [Unref100](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) and [Protein DB](ftp://ftp.wwpdb.org/pub/pdb/derived_data/)) to find homologs. For searching, [MMseqs2](https://github.com/soedinglab/MMseqs2) is utilized which offers a fast and accurate sequence comparison.

The tool has been tested on a conda (v4.10.1) and docker engine (v20.10.7)

## Tool Composition:

- **main.nf** : main script to be called by the user
- **prokka_annotation.nf**  : process for _de novo_ annotation of input fasta
- **mmseqs2.nf**  : process for running additional annotation of input fasta using input query DB
- **update_prokka.nf**  : process for extending prokka annotation by annotatinos found during MMseqs2 process

## Requirements
To use HyPro you only need Nextflow and Conda or Docker/Singularity installed. 

## Script Usage

After installing for example Nextflow and Conda you can either clone this repository to run HyPro or just use Nextflows pull functionality. We also provide a test genome. You can use different Nextflow configuration profiles in order to run the pipeline on different system settings and configurations (e.g. switching from Conda to Docker as the backend software packaging engine):

```bash
nextflow pull hoelzer-lab/hypro
nextflow run hoelzer-lab/hypro -r 0.0.3 -profile local,conda --fasta ~/.nextflow/assets/hoelzer-lab/hypro/test/data/GCF_000471025.2_ASM47102v2_genomic.fna
```

- **-profile local,conda**  : using conda environments for prokka, MMseqs2 and mygene (see configs/conda.config)
- **-profile local,docker** : using docker containers for prokka, MMseqs2 and mygene (see configs/container.config)

When running HyPro multiple times, the tool will look for an existing DB of type **--database** first. If nothing could be found, it will download/build the specified DB and store it in ``nextflow-autodownload-databases/``. Alternatively, you may give HyPro a path to an existing DB created sometime before. Simply hand it over to the [--customdb parameter](#Program-Handling). In this case, do not forget to specify the DB type that you use in **--database**!

```bash
nextflow run hoelzer-lab/hypro -r 0.0.3 -profile local,conda --fasta ~/.nextflow/assets/hoelzer-lab/hypro/test/data/GCF_000471025.2_ASM47102v2_genomic. --database uniprotokb --customdb some/path/to/uniprotkb
```


It is also possible to run HyPro on more than one input genome by using the [--list parameter](#Program-Handling). Instead of passing a single fasta file to **--fasta**, you specify the path to a csv file with two colums per line (sample id, fasta file path) and set the **--list** flag to ``true``.
```bash
nextflow run hoelzer-lab/hypro -r 0.0.3 -profile local,conda --fasta test/input.csv --list true --database uniprotokb --customdb some/path/to/uniprotkb
```

### Program Handling

|Parameter|Type|Description|
|-----|----|-----------|
|**--help** ||Show help message and exit.|
|**--fasta**|String|Path to input genome fasta that shall be annotated. Input can also be a list of multiple fasta files formatted as a .csv file with two columns (sample id and file path)(see **--list**).|
|**--list**|Boolean|Specify whether input is given as a list of files. Default: false|   
|**--database**|String|Specify the target db to search for annotation extension. Current available options: uniprotkb, uniref50, uniref90, uniref100, pdb. Note, that searching on uniref DBs will significantly extend runtime of HyPro. Default: uniprotkb|
|**--custom-db**|String|Specify a path to an existing DB. If no DB is found, HyPro will build it. Requires an according **--database** configuration.|
|**--output**|String|Specify PATH to a directory. HyPro will generate the output structure to PATH. Default: ``results``|
|**--modus**|String|Choose the modus of HyPro to search all hypothetical proteins (full) or leave those out which gained partial annotation (restricted). The dinstinction of fully un-annotated and partial annotated hypothetical proteins was observed for uniprot annotations. Options: full (default), restricted|
|**--threads**|Integer|Define the number of threads to use by MMseqs search and convertalis. Default: 1|
|**--prokka**|String|Control parameters for prokka,e.g. if running HyPro on a bacteria genome that does not follow the standard code.|

### Alignment Parameters

|Parameter|Type|Description|
|----|----|-----------|
|**--evalue**|Float|Include sequence matches with < e-value threshold into the profile. Requires a FLOAT >= 0.0. Default: 0.1|
|**--min-aln-len**|Integer|Specify the minimum alignment length as INT in range 0 to MAX aln length. Default: 0|
|**--pident**|Float|List only matches above this sequence identity for clustering. Enter a FLOAT between 0 and 1.0. Default: 0.0|


## Output

HyPro loads all necessary data for the extension process automatically. It stores all needed information in the **--output** ``PATH``.
For each input fasta it creates a folder with the following files and directories:

* ``prokka.tar.gz`` - stores the prokka annotation for your input genome
* ``mmseqs2_run_*/`` - to store MMseqs2 output and extended prokka annotation. The folder includes one file and two subdirectories:
  * ``mmseqs2_outs/`` storing the ``mmseqs search`` results in tab-separated format (one is the MMseqs2 output while the other contains bit-score-filtered unique hits)
  * ``prokka_restored_updated/`` - all extended files from prokka will be stored here (currently: gff, ffn, faa, gbk)

Note: HyPro will save the MMseqs outputs in blast-like format (tsv) with a unique name composed of the DB you used and the chosen alignment parameters. For example: ``mmseqs2_out_dbuniprotkb_e0.1_a0_p0.0.tsv`` and ``mmseqs2_out_dbuniprotkb_e0.1_a0_p0.0_unique.tsv``will be stored in ``mmseqs2_run_dbuniprotkb_e0.1_a0_p0.0/mmseqs2_outs``. This means the results of an MMseqs2 run on the UniProtKB DB with an e-value cut-off set to 0.1, minimum alignment length of 0 nt and percent identitiy equal to 0 % (those are the default alignment parameters).

A summary file containing the most relevant information on the latest HyPro run with the given MMseqs2 parameterization is also stored in **--output** ``PATH``.

Additionally, the log files for each HyPro process and nextflow execution reports are saved in ``nextflow-run-infos/``. For processes that need to be run with every input sample, the log files are structured into folders named after the input fasta.


## Third-party tools

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|pandas|0.25.2|Might also work for other versions|
|mygene|3.1.0|Automatically installed when running HyPro on a conda or docker engine. Might also work for other versions.|
|mmseqs2|10.6d92c|Automatically installed when running HyPro on a conda or docker engine. <br>**DEPRECATED**: Install in conda environment|
|prokka (recommended)|1.14.6|Used for _de novo_ annotation of test data of chlamydia. Automatically installed when running HyPro on a conda or docker engine.|

