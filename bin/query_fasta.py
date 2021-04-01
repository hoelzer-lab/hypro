#!/usr/bin/env python3

####### PACKAGES ##########
import re, os, sys, subprocess, argparse
# from biothings_client import get_client
import pandas as pd
# from Bio import SeqIO
#######################################################

###### ############GLOBAL VAR #######################
gff_content = {}        # gff content per row number
HyProt_loc = {}         # recognizes the line of each HyProt via prokka feature ID
HyProt_content = {}            # recognizes prokka features of 'hypothetical proteins' via prokka feature ID
################## ARGPARSE #####################
parser = argparse.ArgumentParser()

# Mode
parser.add_argument('-m', '--modus', dest='modus', metavar=['restricted', 'full'], choices=['restricted','full'], default = 'full', required = False, help="Modus of HyPro to decide either for an all hypothetical protein annotation or restricted (only full blanks with no partial annotation). Valid arguments: 'full' and 'restricted'")
#Input-GFF
parser.add_argument('-gff', '--input_gff', dest='input_gff', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the gff file, that shall be extended.")
#Input-FFN
parser.add_argument('-ffn', '--input_ffn', dest='input_ffn', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the ffn file")

args = parser.parse_args()


if isinstance(args.modus, list):
    mode = args.modus[0]
elif isinstance(args.modus, str):
    mode = args.modus

print(f'Start HyPro in {mode} mode')
in_gff = args.input_gff[0]
in_ffn = args.input_ffn[0]

count = 0


################################     FUNCTIONS    #######################################



def load_gff(input=in_gff):
    '''
    Load .gff file from prokka output and read content into global gff_content.
    Check for hypothetical protein entries: store relevant information and keep track of number.
    '''
    row = 0
    HyProt_count = 0
    with open(input, 'r') as gff:
        for line in gff:
            line = line.rstrip('\n')
            attr = line.split('\t')
            gff_content.update({row:attr})      # ggf content with line ID as identifier
            if is_HyProt(attr, 'hypothetical protein'):
                HyProt_count += 1
                save_HyProt(attr, row)
            else:
                pass
            row += 1             # next line
    print("loading gff file finished.")
    print(f"Total Features:\t{row} ")
    print(f"Hypothetical Proteins count:\t{HyProt_count}")
    print(f"Saved records:\t{len(HyProt_content.keys())}")
    return gff_content, HyProt_content   # DEPRECATED

def is_HyProt(attr, regex = ''):
    '''Checks the attribute column for the regex "hypothetical protein" '''
    regex = re.compile(regex)
    try:
        if regex.search(attr[8]):
            return True
        else:
            return False
    except IndexError:
        pass

def save_HyProt(attr, row):
    '''get HyProt IDs, their attributes, and where in gff_content they can be found'''
    global HyProt_loc, HyProt_content, count, mode
    fields = attr[8].split(';')
    ID = fields[0][3:]
    content = []
    # print(f'{ID} {len(fields)} {fields[len(fields)-1]}')
    if mode == 'restricted':
        if len(fields) == 4:                                # hypothetical proteins with no information at all
            for i in range(1, len(fields)):
                content.append(fields[i])
            # count += 1
            HyProt_content.update({ID:content})
            HyProt_loc.update({ID:row})
        else:
            pass
        assert  len(HyProt_content.keys()) == len(HyProt_loc.keys())
    elif mode == 'full':
        if len(fields) == 4 or len(fields) == 5:            # hypothetical proteins with no info and partial info
            for i in range(1, len(fields)):
                content.append(fields[i])
            # count += 1
            HyProt_content.update({ID:content})
            HyProt_loc.update({ID:row})
        else:
            pass
        assert  len(HyProt_content.keys()) == len(HyProt_loc.keys())
    # print(f"save_HyProts:\t{count}")
    return HyProt_content, HyProt_loc   # DEPRECATED

def query_fasta(infile_gff = in_gff, infile_ffn = in_ffn):
    '''
    Match hypothetical protein content from .gff with .ffn file and write information to fasta.
    Resulting file is used by mmseqs2.
    '''
    global HyProt_content, DIR, BN
    print("Try to build query fasta...")
    if bool(HyProt_content):
        fasta = load_fasta(in_ffn)
        with open("query.fasta", 'w') as query:
            for ID in fasta.keys():
                # print(ID)
                if ID in HyProt_content.keys():
                    # print("HyProt")
                    query.write('>' + ID + '\n')
                    query.write(fasta[ID] + '\n')
                else:
                    pass
                    # print('other')
            query.close()
    else:
        print("No query seq saved from input gff. Maybe file is corrupted. Exiting extension pipeline...")
        exit()
    num =  int(subprocess.check_output("grep -c '>' " + "query.fasta", shell=True))
    print(f"Extracted {num} sequences to query.fasta")
    # print(len(HyProt_content.keys()))
    assert num == len(HyProt_content.keys())

def load_fasta(file):
    ''' Load fasta into dict '''
    fasta = {}      # header:seq, all sequences
    header = ''
    seq = ''
    with open(file, 'r') as ffn:
        for line in ffn:
            if line.startswith('>'):
                elem = line.split(' ')
                if header != '' and seq != '':      # update only meaningful seqs
                    fasta.update({header:seq})
                header = elem[0][1:]
                seq = ''
            elif re.match("[ACGT]",line):           # nt seq of hypothetical protein
                seq = seq + str(line).rstrip('\n')
    fasta.update({header:seq})
    # print(len(fasta.keys()))
    return fasta
    # print(fasta)


############# MAIN ######################

load_gff()
query_fasta()
