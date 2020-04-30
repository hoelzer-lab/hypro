<!--[![Build Status](https://travis-ci.org/hoelzer/ribap.svg?branch=master)](https://travis-ci.org/hoelzer/ribap)-->
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Python](https://img.shields.io/badge/Language-Python3.7-green.svg)
[![Twitter Follow](https://img.shields.io/twitter/follow/martinhoelzer.svg?style=social)](https://twitter.com/martinhoelzer) 

Authors: Martin H&ouml;lzer, Maximilian Arlt

# HyPro
Protein-coding annotation extension using additional homology searches against larger databases.

## Summary

The HyPro tool extends common protein-coding annotations made with [Prokka](https://github.com/tseemann/prokka) using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (available are [UniProtKB](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/), [Uniref50](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/), [Uniref90](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/), [Unref100](ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/) and [Protein DB](ftp://ftp.wwpdb.org/pub/pdb/derived_data/)) to find homologs. For searching, [MMseqs2](https://github.com/soedinglab/MMseqs2) is utilized which offers a fast and accurate sequence comparison.

The tool has been tested in a conda environment (v4.7.11).

## Tool Composition:

- **hypro.py**     main script to be called by the user
- **mmseqs2.sh**     bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting

## Requirements
HyPro requires the provided list of software to function properly. 

|Program/Package|Version|Note|
|---------------|-------|------|
|python|3.7|Might also work for other python3 versions|
|mygene|3.1.0|Might also work for other versions|
|pandas|0.25.2|Might also work for other versions.|
|mmseqs2|10.6d92c|Install in conda environment|
|prokka (recommended)|1.14.6|used for _de novo_ annotation of test data of chlamydia|

It is recommended to clone this repository and use a conda environment for HyPro.


```bash
git clone https://github.com/hoelzer-lab/hypro.git
cd hypro
```
Create a conda environment:

```
conda create -n hypro python=3.7 pandas=0.25.2 mmseqs2=10.6d92c prokka=1.14.6 mygene=3.1.0
conda activate hypro
```

Or simply install HyPro from the bioconda repository:

```
conda create -n hypro python=3.7
conda activate hypro
conda install -c bioconda hypro
```

## Script Usage

After installing all dependencies (see commands above) and cloning this repository we simply use Prokka on a test genome: 

```bash
prokka --prefix testrun --outdir run/prokka test/data/GCF_000471025.2_ASM47102v2_genomic.fna
```

and then execute the scripts from this repository:

```bash
scripts/hypro.py -i run/prokka/testrun.gff -o run/hypro -d uniprotkb -f scripts/mmseqs2.sh -m full
```
assuming that your current working directory is the previously downloaded git repository of HyPro. Otherwise, please adjust ``scripts/`` accordingly.
or in case of conda package simply do:
```bash
hypro.py -i run/prokka/testrun.gff -o run/hypro -d uniprotkb -f mmseqs2.sh -m full
```
No path to the ``mmseqs2.sh`` is required.

When running HyPro to the same output folder multiple times, the tool will look for an existing DB of the given type ('-d') first. If nothing could be found, it will download/build the specified DB. Alternatively, you may give HyPro a path to an existing DB created sometime before. Simply hand it over to the [-c parameter](#Program-Handling). In this case, do not forget to specify the DB type that you use in -d!

```bash
scripts/hypro.py -i run/prokka/testrun.gff -o run/hypro_re-use_db -f scripts/mmseqs2.sh -m full -d uniref50 -c run/hypro/db/uniref50
```

### Program Handling
|Short|Long|Description|
|-----|----|-----------|
|**-h**|**--help** |Show this help message and exit.|
|**-i**|**--input**|Path to input gff that shall be extended.|   
|**-o**|**--output**|Specify PATH to a directory. HyPro will generate the output structure to PATH.|
|**-d**|**--database**|Specify the target db to search for annotation extension. Current available options: uniprotkb, uniref50, uniref90, uniref100, pdb. Note, that searching on uniref DBs will significantly extend runtime of HyPro. [uniprotkb]|
|**-f**|**--mmseq2**|Specify the path to the mmseqs2.sh. If using the conda package, 'mmseqs2.sh' is enough.|
|**-m**|**--modus**|Choose the modus of HyPro to search all hypothetical proteins (full) or leave those out which gained partial annotation (restricted). The dinstinction of fully un-annotated and partial annotated hypothetical proteins was observed for uniprot annotations. Options: [full, restricted]|
|**-c**|**--custom-db**|Specify a path to an existing DB. If no DB is found, HyPro will build it. Requires an according -d configuration.|
|**-t**|**--threads**|Define the number of threads to use by mmseqs indexdb, search and convertalis. [1]|

### Alignment Parameters
|Short|Long|Description|
|-----|----|-----------|
|**-e**|**--evalue**|Include sequence matches with < e-value threshold into the profile. Requires a FLOAT >= 0.0. [0.1]|
|**-a**|**--min-aln-len**| Specify the minimum alignment length as INT in range 0 to MAX aln length. [0]|
|**-p**|**-pident**|List only matches above this sequence identity for clustering. Enter a FLOAT between 0 and 1.0. [0.0]|

## Output

HyPro loads all necessary data for the extension process automatically. It stores all needed information in the ``-o/--output PATH``. In ``PATH``, it creates the following directories:

* ``db`` - stores the databases you have chosen, each in an own directory
* ``mmseqs2_output`` - to store mmseqs2 output. The folder includes a subdirectory "final_outs" storing the ``mmseqs search`` results in tab-separated format (one is the mmseqs2 output while the other contains bit-score-filtered unique hits)
* ``output`` - all extended files from prokka will be stored here (currently: gff, ffn, faa, gbk)

Note: HyPro will save the mmseqs outputs in blast-like format (tsv) with a unique name composed of the DB you used and the chosen alignment parameters. For example: *mmseqs2_out_db_uniprotkb_e0.1_a0_p0.0.tsv* and *mmseqs2_out_db_uniprotkb_e0.1_a0_p0.0_unique.tsv*. This means the results of an mmseqs run on the UniProtKB DB with an e-value cut-off set to 0.1, minimum alignment length of 0 nt and percent identitiy equal to 0 % (those are the default alignment parameters).
