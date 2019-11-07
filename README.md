## prokkaX
Protein-coding annotation extension using additional homology searches against larger databases.

# Summary

The prokkaX tool extends common protein-coding anotations made with Prokka using additional homology searches. The approach currently takes a gff input file, extracts the sequences of hypothetical proteins and searches against a selected database (currently available: only uniprotkb) to find homologs. For searching, mmseqs2 is utilized which currently offers a fast as well as accurate sequence comparison, to the best of our knowledge.

The tool has been tested in an anaconda environment (v. 4.7.11). 

-----------
# Tool Composition:

*gff_extend.py* -  _main script - to be called by the user_
*mseqs2.sh*     - _bash script comprising all required mmseqs2 commands (createdb, createindex, search) and output formatting_

-----------

# Requirements
ProkkaX requires the provided list of software to function properly. It is recommended to use a conda environment for prokkaX.

|Program/Package|Version|Note|
|---------------|-------||
|python|3.7|Might also work for other python3 versions|
|pandas|||
|mmseqs2|10.6d92c||
||||
gff_content.py -i PATH -o PATH -d STR -m PATH

optional arguments:
  -h, --help              show this help message and exit
  
  -i PATH, --input PATH   Path to input gff that shall be extended
  
  -o PATH, --output PATH  Path to an directory. The program will build all output directories in here.
  
  -d STR, --database STR  Specifiy target db to use for extension. Available
                          options are 'uniprotkb'.
                          
  -m PATH, --mmseq2 PATH  Specify the path to the mmseqs2.sh.
