# prokkaX
Extend common protein-coding annotations made w/ Prokka using additional homology searches against larger databases

# Main script

gff_content.py -i PATH -o PATH -d STR -m PATH

optional arguments:
  -h, --help              show this help message and exit
  
  -i PATH, --input PATH   Path to input gff that shall be extended
  
  -o PATH, --output PATH  Path to an directory. The program will build all output directories in here.
  
  -d STR, --database STR  Specifiy target db to use for extension. Available
                          options are 'uniprotkb'.
                          
  -m PATH, --mmseq2 PATH  Specify the path to the mmseqs2.sh.
