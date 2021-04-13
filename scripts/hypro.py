#!/usr/bin/env python3

####### PACKAGES ##########
import re, os, sys, subprocess, argparse
import mygene
# from biothings_client import get_client
import pandas as pd
# from Bio import SeqIO
#######################################################

###### ############GLOBAL VAR #######################
gff_content = {}        # gff content per row number
HyProt_loc = {}         # recognizes the line of each HyProt via prokka feature ID
HyProt_content = {}            # recognizes prokka features of 'hypothetical proteins' via prokka feature ID
valid_db = ['uniprotkb','uniref100', 'uniref90', 'uniref50', 'pdb']        # to be extended for novel db-types
id_scinames = {}               # feature IDs and target scientific name
automated_dbload = True
################## ARGPARSE #####################
parser = argparse.ArgumentParser()

# Mode
parser.add_argument('-m', '--modus', dest='modus', metavar=['restricted', 'full'], choices=['restricted','full'], default = 'full', required = False, help="Modus of HyPro to decide either for an all hypothetical protein annotation or restricted (only full blanks with no partial annotation). Valid arguments: 'full' and 'restricted'")
#Input-GFF
parser.add_argument('-i', '--input', dest='input', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the gff file, that shall be extended.")
# Output_dir
parser.add_argument('-o', '--output', dest='output', action='store', metavar='PATH', nargs=1, required=True, default='..', help="Specify PATH to a directory. HyPro will generate all output to this.")
# DB-Type
parser.add_argument('-d', '--database', dest='db', action='store', metavar='STR', nargs=1, required=False, default='uniprotkb', help="Specify the target db to search in for annotation extension. Available options: 'uniprotkb', 'uniref100', 'uniref90', 'uniref50', 'pdb' [uniprotkb]")
# location of mmseq2.sh:
parser.add_argument('-f', '--mmseqs2', dest='mmseq', action='store', metavar='PATH', nargs=1, required=False, default="./scripts/mmseqs2.sh", help="Specify the path to the mmseqs2.sh. If using the conda package, 'mmseqs2.sh' is enough. Default path is './scripts/msmeqs2.sh'")
# ALL or just BLANKS ?
#parser.add_argument('-m', '--modus', dest='modus', action='store', metavar='STR', nargs=1, default='blanks', help="Search all hypothetical proteins or just the blanks. default: blanks")
# Custom-DB
parser.add_argument('-c', '--custom-db', dest='custdb', action='store', metavar='STR', nargs=1, required=False, help="Specify a path. If no database is found, HyPro will build it.  Requires an according -d configuration.")
# E-Value
parser.add_argument('-e', '--evalue', dest='evalue', action='store', metavar='FLOAT', nargs=1, required=False, default='0.1', help='Include sequence matches with < e-value threshold into the profile. Requires a FLOAT >= 0.0. [0.1]')
# Alignment length
parser.add_argument('-a', '--min-aln-len', dest='alnlen', action='store', metavar='INT', nargs=1, required=False, default='0', help='Specify the minimum alignment length as INT in range 0 to MAX aln length. [0]')
# Percentage identity
parser.add_argument('-p', '--pident', dest='pident', action='store', metavar='FLOAT', nargs=1, required=False, default='0.0', help='List only matches above this sequence identity for clustering. Enter a FLOAT between 0 and 1.0. [0.0]')
# sensitivity of mmseqs2
# parser.add_argument('-s', '--sens', dest='sens', action='store', metavar='FLOAT', nargs=1, required=False, default='5.7', help=' Adjust the sensitivity of mmseqs search/index. 1.0 fastest, 7.5 most sensitive (range: 1.0 - 7.5). [5.7]')
# Threads
parser.add_argument('-t', '--threads', dest='threads', action='store', metavar='INT', nargs=1, required=False, default='1', help='Define number of threads to use by mmseqs indexdb, mmseqs search and mmseqs convertalis. [1]')

args = parser.parse_args()

# if isinstance(args.sens, list):
#     sens = float(args.sens[0])
# elif isinstance(args.sens, str):
#     sens = float(args.sens)

if isinstance(args.mmseq, list):
    ms = args.mmseq[0]
elif isinstance(args.mmseq, str):
    ms = args.mmseq

if isinstance(args.threads, list):
    threads = int(args.threads[0])
elif isinstance(args.threads, str):
    threads = int(args.threads)

if bool(args.custdb):               # set
    automated_dbload = False
    custdb = args.custdb[0]

if isinstance(args.modus, list):
    mode = args.modus[0]
elif isinstance(args.modus, str):
    mode = args.modus

print(f'Start HyPro in {mode} mode')
in_gff = args.input[0]
out_dir = args.output[0]

if isinstance(args.db, list):
    db = args.db[0]
elif isinstance(args.db, str):
    db = args.db

print(out_dir)

# ms = args.mmseq[0]

if isinstance(args.evalue, list):
    evalue=float(args.evalue[0])
elif isinstance(args.evalue, str):
    evalue=float(args.evalue)

if isinstance(args.alnlen, list):
    alnlen=int(args.alnlen[0])
elif isinstance(args.alnlen, str):
    alnlen=int(args.alnlen)

if isinstance(args.pident, list):
    pident=float(args.pident[0])
elif isinstance(args.pident, str):
    pident=float(args.pident)

DIR = ''            # input directory, all prokka files should reside at
BN = ''             # basename of prokka output files
EXT = ''            # extension found in input
count = 0

################# FUNCTIONS ########################


def update_faa(output): #former 2nd arg: in_gff
    global DIR, BN, id_scinames
    faa_list = []
    infile = f'{DIR}/{BN}.faa'  # input faa
    # print(id_scinames)
    with open(infile, 'r') as faa:
        for line in faa:
            if line.startswith('>'):
                elem = line.split(' ')
                if elem[0][1:] in id_scinames.keys():
                    # print(elem[0][1:])
                    # print(id_scinames[elem[0][1:]])
                    header = f'{elem[0]} {id_scinames[elem[0][1:]]}\n'          # new header
                    faa_list.append(header)                                     # update header
                    # print(f'{elem[0]} \t {line} \t {header}')
                else:
                    faa_list.append(line)
            else:
                faa_list.append(line)

    #write the updated faa file
    with open(output + '/output/' + BN + "_extended.faa", 'w') as faa:
        for entry in faa_list:
            faa.write(entry)

def update_ffn(output): #former 2nd arg: in_gff
    global DIR, BN, id_scinames
    ffn_list = []
    infile = f'{DIR}/{BN}.ffn'  # input faa
    # print(id_scinames)
    with open(infile, 'r') as ffn:
        for line in ffn:
            if line.startswith('>'):
                elem = line.split(' ')
                if elem[0][1:] in id_scinames.keys():
                    # print(elem[0][1:])
                    # print(id_scinames[elem[0][1:]])
                    header = f'{elem[0]} {id_scinames[elem[0][1:]]}\n'          # new header
                    ffn_list.append(header)                                     # update header
                    # print(f'{elem[0]} \t {line} \t {header}')
                else:
                    ffn_list.append(line)
            else:
                ffn_list.append(line)
    # print(faa_list)

    #write the updated faa file
    with open(output + '/output/' + BN + "_extended.ffn", 'w') as faa:
        for entry in ffn_list:
            faa.write(entry)

# extend gbk
def extend_gbk(output, id_alninfo):
    '''extend the information in prokka output gbk with mmseq2 homology findings - based on the description of genbank file structure of NCBI'''
    global DIR, BN, id_scinames
    infile = f'{DIR}/{BN}.gbk'          #create the genbank file name from loaded DIR path and file basename
    updated = 0                         # assert to check all HyProts had been matched
    with open (infile, 'r') as file:
        line = "\n"                     # initialize non-empty string to pass the first while statement
        content = []
        while bool(line):               # as long as line is not empty = reading file is not finished
            comment = []                            # comment section content
            features = []                           # the feature section content
            end_info = []                           # the ORIGIN part; after feature, not changed
            while "COMMENT" not in line:
                content.append(line)
                # ahead_info.append(line)
                line = file.readline()
            # print(line)
            while "FEATURES" not in line:
                comment.append(line)
                line = file.readline()
            comment = update_comment(comment)
            content.extend(comment)
            while "ORIGIN" not in line:
                if line[5] != ' ':              # acc. to def. of genbank file format in NCBI: feature starts in 6th column
                    # feature = []          # one complete feature found by prokka
                    global lspaces          # global to use in format_notes(), too
                    descriptors = {}        # key: descriptor, value: content (e.g. "/product":"hypothetical protein")
                    translation_ext = []        # from 2nd line of translation descriptor, until the end: first line parsed to dict
                    features.append(line)
                    line = file.readline()
                    lclean_line = line.lstrip()         # clip the leading whitespaces of the line
                    lspaces = len(line) - len(lclean_line)    # get the number of white spaces
                    while line[5] == ' ':               # as long as descriptors of the same feature parsed
                        # print(line)
                        if lclean_line.startswith('/'):
                            elem = lclean_line.split('=')       # split the feature content line
                            if elem[0] in descriptors:
                                descriptors[elem[0]].append(f'{elem[0]}={elem[1]}')
                            else:
                                descriptors[elem[0]] = [f'{elem[0]}={elem[1]}']
                        else:                                               # only 2nd-to-last translation lines shall be parsed in here
                            descriptors[elem[0]].append(lclean_line)        # when feature is complete, check for update in mmseq hit dictionary; update if found
                        line = file.readline()
                        lclean_line = line.lstrip(' ')  # clip whitespaces
                # Add info of the feature, if it was build:
                if bool(descriptors):
                    try:
                        ID = descriptors['/locus_tag'][0][11:]
                        ID = ID.strip()
                        ID = ID.strip('\"')
                        if ID in id_scinames.keys():       # look for mmseq2 extension
                            # print("Updating a match!")
                            updated += 1                    # count number of updated features
                            descriptors = insert_mmseq_info(descriptors, id_alninfo)
                            id_scinames.pop(ID)
                    except KeyError:                    # if no descriptor '/locus tag' was found, the loaded info is not a real feature (source or anythin else...)
                        pass
                    # print(descriptors)
                    features.append(descriptors)
                if bool(translation_ext):
                    features.append(translation_ext)

            # assert len(id_scinames.keys()) == updated
            content.extend(features)
            while "//" not in line:
                end_info.append(line)
                line = file.readline()
            end_info.append(line)
            content.extend(end_info)
            line = file.readline()

        write_gbk(output, content, lspaces)
        print(f'Updated {updated} features in the genbank file.')

def format_notes(string, delimiter):
    '''Format a to lines of 80 characters which corresponds to gbk convention'''
    global lspaces
    elem = string.split(delimiter)              # split the string into its meaningful units
    l_string = ""
    lines = 0
    # len_del = len(delimiter)                  # if sth else than a single sign is used to delimit
    for qual in elem:
        if lines == 0:                          # add the qualifier after the first leading spaces
            l_string += '/notes=\"'
        else:
            l_string += (lspaces * ' ')
        lines += 1
        l_string += qual + '\n'
    l_string = l_string.rstrip('\n')
    l_string += '\"\n'
    return l_string

def insert_mmseq_info(descriptors, id_alninfo):
    global id_scinames, db
    ID = descriptors['/locus_tag'][0].strip()
    ID = ID[11:]
    ID = ID.strip('\"')
    descriptors['/product'][0] = f'/product="{id_scinames[ID]}"\n'
    descriptors['/gene'] = [f'/gene="{id_alninfo[ID][1]}"\n']
    descriptors['/inference'].append(f'/inference="mmseqs2 {db}"\n')
    alninfo = id_alninfo[ID][0].replace("=",':')                    # change the separator from "=" (used in gff) to ":", since "=" is the separator for the descriptor lines
    alninfo = alninfo.rstrip(';')
    alninfo = format_notes(alninfo, ';')

    descriptors['/notes'] = [alninfo]               # save the alignment info with separator ':'
    # print(descriptors['/notes'])
    return descriptors

def update_comment(in_list):
    add_string = ''
    w_spaces = len(in_list[1]) - len(in_list[1].lstrip(' '))
    add_string = (w_spaces * " ") + "extended by homology searches with HyPro v0.1 (based on mmseq2)\n"
    in_list.append(add_string)
    return in_list

def write_gbk(output, content, wspaces):
    global DIR, BN
    with open(output + "/output/" + BN + "_extended.gbk", 'w') as gbk_out:
        for elem in content[1:]:
            # if elem == content[0]:
            #     pass
            if isinstance(elem, dict):  # write the updated features
                for key in elem.keys():
                    if "translation" not in key:                                # write translation the last
                        for val in elem[key]:
                            gbk_out.write(wspaces * ' ' + str(val))             # write the key array values one by one
                    else:
                        pass
                try:
                    for val in elem['/translation']:
                        gbk_out.write(wspaces * ' ' + str(val))   # ensures to write the first line of the translated protein the last!
                except KeyError:
                    pass
            elif isinstance(elem, list):
                for place in elem:
                    gbk_out.write(place)
            else:
                gbk_out.write(elem)


############ update_gff
def rusure(dbtype):
    if dbtype == valid_db[3]:
        print("Downloading the uniref50 may take several hours. Are you sure to continue? [yes/no]")
    elif dbtype == valid_db[2]:
        print("Downloading the uniref90 might take several days. Are you sure to continue? [yes/no]")
    elif dbtype == valid_db[1]:
        print("Downloading the uniref100 might take several days. Are you sure to continue? [yes/no]")

    x = input()
    if x == 'yes':
        pass
    elif x == 'no':
        print('You can choose another database when re-running. Exiting...')
        exit()
    else:
        print('Unknown option specified. Please choose: [yes/no]')
        rusure()

def load_gff(input=in_gff):
    '''
    Load .gff file from prokka output and read content into global gff_content.
    Check for hypothetical protein entries: store relevant information and keep track of number.
    '''
    # is_outdir()
    # is_gff()
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
    # print(f"HyProt_content:\t{len(HyProt_content.keys())}")
    # print(gff_content[len(gff_content)-1])
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

def query_fasta(infile = in_gff, output = out_dir):
    '''
    Match hypothetical protein content from .gff with .ffn file and write information to fasta.
    Resulting file is used by mmseqs2.
    '''
    global HyProt_content, DIR, BN
    # print("Build output directory 'HyPro' in specified output location, if not existing.")
    print("Try to build query fasta...")
    ffn_file = DIR + '/' + BN + '.ffn'      # feature sequences
    if bool(HyProt_content):
        fasta = load_fasta(ffn_file)
        with open(output + "/query.fasta", 'w') as query:
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
    num =  int(subprocess.check_output("grep -c '>' " + output +  "/query.fasta", shell=True))
    print(f"Extracted {num} sequences to {output}/query.fasta")
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

# download a particular db
def download_db(output = out_dir, dbtype = db): #--------------------------TO BE EXTENDED WITH ADDITIONAL DBs-----------------------------------
    global valid_db, automated_dbload
    weblink = ''
    if automated_dbload:                        # if no custom path defined, define db path as the output path
        path = output.rstrip('/') + "/db/" + dbtype
    else:
        global custdb
        path = custdb
    if dbtype == valid_db[0]:                   # if clause to decide which db to download
        db_fasta = path + "/uniprot_sprot.fasta"
        db_target = path + "/target_db"
        if loaded(db_fasta, automated_dbload):
            print("Uniprotkb fasta exists already. Skip download.")
        else:
            print("No db found. Download...")
            weblink = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
            file = "uniprot_sprot.fasta.gz"
            print(path)
            pwd = os.getcwd()
            # print(pwd)
            os.chdir(path)
            # print(os.getcwd())
            print(f"Download database to {path}")
            os.system(f"wget {weblink}")
            os.system(f"gunzip ./{file}")
            os.chdir(pwd)
            # print(os.getcwd)
    elif dbtype == valid_db[1]:
        db_fasta = path + "/uniref100.fasta"
        db_target = path + "/target_db"
        if loaded(db_fasta, automated_dbload):
            print("Uniref100 fasta exists already. Skip download.")
        else:
            rusure(dbtype)
            print("No db found. Download...")
            weblink = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
            file = "uniref100.fasta.gz"
            print(path)
            pwd = os.getcwd()
            # print(pwd)
            os.chdir(path)
            # print(os.getcwd())
            print(f"Download database to {path}")
            os.system(f"wget {weblink}")
            os.system(f"gunzip ./{file}")
            os.chdir(pwd)
            # print(os.getcwd)
    elif dbtype == valid_db[2]:
        db_fasta = path + "/uniref90.fasta"
        db_target = path + "/target_db"
        if loaded(db_fasta, automated_dbload):
            print("Uniref90 fasta exists already. Skip download.")
        else:
            rusure(dbtype)
            print("No db found. Download...")
            weblink = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
            file = "uniref90.fasta.gz"
            print(path)
            pwd = os.getcwd()
            # print(pwd)
            os.chdir(path)
            # print(os.getcwd())
            print(f"Download database to {path}")
            os.system(f"wget {weblink}")
            os.system(f"gunzip ./{file}")
            os.chdir(pwd)
            # print(os.getcwd)
    elif dbtype == valid_db[3]:
        db_fasta = path + "/uniref50.fasta"
        db_target = path + "/target_db"
        if loaded(db_fasta, automated_dbload):
            print("Uniref50 fasta exists already. Skip download.")
        else:
            rusure(dbtype)
            print("No db found. Download...")
            weblink = "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
            file = "uniref50.fasta.gz"
            print(path)
            pwd = os.getcwd()
            # print(pwd)
            os.chdir(path)
            # print(os.getcwd())
            print(f"Download database to {path}")
            os.system(f"wget {weblink}")
            os.system(f"gunzip ./{file}")
            os.chdir(pwd)
            # print(os.getcwd)
    elif dbtype == valid_db[4]:
        db_fasta = path + "/pdb_seqres.txt"
        db_target = path + "/target_db"
        if loaded(db_fasta, automated_dbload):
            print("Uniref50 fasta exists already. Skip download.")
        else:
            rusure(dbtype)
            print("No db found. Download...")
            weblink = "ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz"
            file = "pdb_seqres.txt.gz"
            print(path)
            pwd = os.getcwd()
            # print(pwd)
            os.chdir(path)
            # print(os.getcwd())
            print(f"Download database to {path}")
            os.system(f"wget {weblink}")
            os.system(f"gunzip ./{file}")
            os.chdir(pwd)
            # print(os.getcwd)
    else:
        print("Error. Exit from download_db function...")
        exit()
    return path, db_fasta, db_target    #-----------------------------------------------------------------------------------------

def mmseq(dbfasta, dbtarget, output=out_dir, mmseq=ms, ev=evalue, aln=alnlen, pi=pident, t=threads):
    '''
    Search hypothetical proteins from prokka output in a selected database using mmseqs2.
    Extract, format and write new annotation to output id_infos, update previous HyProt_content.
    '''
    global HyProt_content, db, id_scinames
    id_infos = {}                                           # additional information found  with mmseq; saved with ID as key
    output = output.rstrip('/') + '/mmseqs_output'
    # # execute mmseq2
    os.system(f"{mmseq} {output} {dbfasta} {dbtarget} {db} {ev} {alnlen} {pident} {t}")

    # # fraction identified
    hit_nums = str(subprocess.check_output(f"cut -f1 {output}/final_outs/mmseqs2_out_db_{db}_e{ev}_a{alnlen}_p{pident}_unique.tsv | wc -l", shell=True))
    print(f"The homology search assigns a function to {int(hit_nums[2:-3])-1} / {len(HyProt_content.keys())} hypothetical proteins\n\t(this includes also hits on hypothetical proteins in the database)")

    # load annotation pandas dataframe
    mmseqs_out = pd.read_csv(f"{output}/final_outs/mmseqs2_out_db_{db}_e{ev}_a{alnlen}_p{pident}_unique.tsv" , sep='\t')

    #create dict from pandas df - ID:annotation dictionary
    out_dict = mmseqs_out.to_dict('index')  # dict of dicts; {index:{row}}; row = {col-X:row-entry}
    # if out_dict is unref50/90/100:
    if db == valid_db[1] or db == valid_db[2] or db == valid_db[3]:
        for num in out_dict.keys():
            # out_dict[num]["target"]
            value = out_dict[num]
            target_split = out_dict[num]['target'].split('_')
            out_dict[num].update({'target':target_split[1]})  # ID form "Uniref[50,90,100]_uniprotID" replaced by uniprotID
            # print(out_dict[num]['target'])

    # print(f'MMSeqs2 Output:\n{out_dict}')
    # print(prots_dict.keys())

    # lookup sci_names and symbols of target IDs
    dict_scinames = get_names(f"{output}/final_outs/mmseqs2_out_db_{db}_e{ev}_a{alnlen}_p{pident}_unique.tsv", db)

    #print(f'SciNames of TargetIDs:\n{dict_scinames}')
    # save findings in HyProt_content (target info + scinames/symbols of target IDs)
    # print(HyProt_content)
    for index in out_dict.keys():
        # print(out_dict[index]['query'])
        HyProt = out_dict[index]['query']
        db_ID = out_dict[index]['target']
        # print(dict_scinames[db_ID][1])
        id_scinames.update({HyProt:dict_scinames[db_ID][1]})                # save scinames under the HyProt prokka IDs
        if HyProt in HyProt_content.keys():
            info = ''
            for attr in out_dict[index].keys():
                info = f"{info}HyPro_{attr}={out_dict[index][attr]};"
            info = f"{info}HyPro_symbol={dict_scinames[db_ID][0]};"      # add target symbol to info
            info = f"{info}HyPro_sciname={dict_scinames[db_ID][1]};"     # add target sci name to info
            id_infos[HyProt] = [info]
            id_infos[HyProt].append(out_dict[index]['target'])
                # print(info)
                # print(out_dict[index][attr])
            # print(info)
            # print(HyProt_content[HyProt])
            HyProt_content[HyProt][0] = f"inference=mmseqs2 {db}"
            HyProt_content[HyProt][2] = f"product={out_dict[index]['target']}"
            HyProt_content[HyProt].append(info.rstrip(';'))
    counts = count_HyProts(id_scinames)
    print(f"Of those: {counts} / {int(hit_nums[2:-3])-1} are hypothetical proteins.")
    return id_infos

def update_gff(output=out_dir, delimiter = '\t'):
    global HyProt_content, HyProt_loc, gff_content, BN
    for HyProt in HyProt_content.keys():
        data = 'ID='
        data = f'{data}{HyProt}'
        # print(HyProt)
        # print(HyProt_content[HyProt])
        for elem in HyProt_content[HyProt]:
            data = f"{data};{elem}"
            # print(data)
        # data = data + '\''
        # print(data)
        loc = HyProt_loc[HyProt]
        # print(gff_content[loc])
        gff_content[loc][8] = data
        # print(gff_content[loc])
    # write to gff
    with open(output + "/" + BN + "_extended.gff", 'w') as gff_up:
        print(f"Write extended gff:\t{output}/{BN}_extended.gff")
        for key in gff_content:
            line = ''
            # print(gff_content[key])
            for elem in gff_content[key]:
                line = line + str(elem) + '\t'
                # print(line)
            # break
            #     print(line)
            gff_up.write(str(line.rstrip('\t')) + '\n')

def get_input_info(input = in_gff):
    ''' Store information about specified input '''
    DIR = os.path.dirname(input)
    NAME = os.path.basename(input)
    file_split = NAME.split('.')
    BN = file_split[0]
    EXT = file_split[1]
    return DIR, BN, EXT

def check_args(infile = in_gff, output = out_dir, mmseq = ms, dbtype = db):
    '''Control all given paths and files first, before running anything!'''
    global DIR, BN, EXT, automated_dbload
    print("CHECK ARGUMENTS")
    print(f"Current working directory: {os.getcwd()}")
    # check input gff file is valid - does not really test gff format!
    DIR, BN, EXT = get_input_info(infile)
    mmseq_name = os.path.basename(mmseq)
    if os.path.isfile(infile) and EXT == 'gff':
        print("Specified input is a file with 'gff' extension.")
    else:
        print("Invalid input file. Please enter '--help' flag for information on script usage.")
        exit()
    # check output directory path
    print(f'Check args: outdir is: {output}')
    print("Building the output structure...")
    create_outdir(output, dbtype)
    print("Done.")

    if dbtype in valid_db:                               # db is allowed entry
        print(f"Extension will be performed on {db}")
    else:
        print("Unknown dbtype chosen. please check '--help' to choose a valid database for mmseqs2.")
        exit()
     # check, if ffn file is present
    ffn_file = DIR + '/' + BN + '.ffn'      # feature sequences
    if os.path.isfile(ffn_file):
        print("Feature Fasta file found.")
    else:
        print(f"{filename} does not exist. To generate mmseqs2 input query, a 'ffn'-file is required. Make sure, your prokka output is not corrupted.")
        exit()

    # check, if mmseq2.sh is valid
    try:
        path_to_sh = str(subprocess.check_output(f'which {mmseq}', shell=True)).rstrip('\n')
    except subprocess.CalledProcessError:
        print("mmseqs2.sh path specification is corrupted. Check the specified path and the file you referenced.")
        exit()

    if os.path.isfile(path_to_sh[2:-3]) and mmseq_name == "mmseqs2.sh":     # path should point to mmseqs2.sh
        print("Found mmseqs2.sh.")
    else:
        print("mmseqs2.sh path specification is corrupted. Check the specified path and the file you referenced.")
        exit()

    # check, if customdb is a valid path
    if automated_dbload == False:
        global custdb
        if os.path.isdir(custdb):
            print("Path to custom db is valid.")
        else:
            print("the defined path to a custom db is corrupted. Please check. Exiting...")
            exit()

def loaded(dbfasta, automated_dbload):
    global valid_db
    if automated_dbload:
        if os.path.isfile(dbfasta):
            return True
        else:                       # Continue script, if automatic download does not find a db fasta
            return False
    else:
        if os.path.isfile(dbfasta):
            return True
        else:                           # exit, if custom db is not find in path. Maybe the db is corrupt or the -d does not match the db in -c. This limits the usage of -c to only customdbs
            print('No valid fasta file of specified DB Type in custom path. Did you specify -d (dbtype) according to -c (path to custom db)?\n Exiting...')
            exit()

def create_outdir(out_dir, dbtype):
    '''create the output path structure'''
    os.system(f'mkdir -p {out_dir}')
    os.system(f'mkdir -p {out_dir}/db')
    os.system(f'mkdir -p {out_dir}/db/{dbtype}')
    os.system(f'mkdir -p {out_dir}/mmseqs_output/tmp')
    os.system(f'mkdir -p {out_dir}/output')

def get_names(table, db):   # table should be an input
    global valid_db
    gene_ids = []
    translate = {}
    ids = []
    with open(table,'r') as csv:
        csv.readline()
        for line in csv:
            elem = line.split('\t')
            gene = mygene.MyGeneInfo()
            ids.append(elem[1])
            # print(gene.getgene('uniport:P24941','name,symbol'))
    # print(f'{len(ids)} {sorted(ids)}')
    #print(ids)
    if db == valid_db[0]:
        ret = gene.querymany(ids, scopes='uniprot', fields='name,symbol', verbose= False)           # set verbose 'True' to get info about duplicates/missing values of name parsing
        sciname_dict = collect_scinames(ret)
    elif db == valid_db[1] or db == valid_db[2] or db == valid_db[3]:
        for id in range(len(ids)):
            entry = ids[id].split('_')
            ids[id] = entry[1]
       # print(ids)
        ret = gene.querymany(ids, scopes='uniprot', fields='name,symbol', verbose= False)           # set verbose 'True' to get info about duplicates/missing values of name parsing
        sciname_dict = collect_scinames(ret)
    elif db == valid_db[4]:
        ret = gene.querymany(ids, scopes='pdb', fields='name,symbol', verbose= False)           # set verbose 'True' to get info about duplicates/missing values of name parsing
        sciname_dict = collect_scinames(ret)
    else:
        print('Sth went wrong in function get_names.')
        exit()

    # print(sciname_dict)
    return sciname_dict
    # ret['query','name','symbol']
    # ret.to_csv('/mnt/mahlzeitlocal/projects/ma_neander_assembly/hiwi/HyPro/test/HyPro/mmseq_output/mmseq2_out_unique_names.tsv')
    # print(ret)
def collect_scinames(name_list):
    sciname_dict = {}               # {ID:[symbol, sciname]}
    for entry in name_list:
        try:
            sciname_dict.update({entry['query']:[entry['symbol'],entry['name']]})
        except:
            sciname_dict.update({entry['query']:['HyProt','hypothetical protein']})
        # print(entry['query'])
    return sciname_dict

def count_HyProts(features):
    HyProt_count = 0
    HyProt = re.compile("hypothetical protein")
    for id in features.keys():
        if re.search(HyProt,features[id]):
            HyProt_count += 1
    return HyProt_count



############# MAIN ######################

check_args()
create_outdir(out_dir, db)
load_gff()
query_fasta()
db_dir, dbfasta, dbtarget = download_db()
id_alninfo = mmseq(dbfasta, dbtarget)
update_gff()
update_faa(out_dir)
update_ffn(out_dir)
extend_gbk(out_dir, id_alninfo)
# get_names('/mnt/mahlzeitlocal/projects/ma_neander_assembly/hiwi/HyPro/test/HyPro/mmseq_output/mmseq2_out_unique.tsv')


# TO DO:
    # object oriented - > atrributes: ID eC_number Name db_xref gene inference locus_tag product
    # design for different db
    # mmseq() - parametric input of target db info


# Requirements

# prokka 1.14.0
# mmseqs2 10.6d92c
# mygene module 3.1.0
# python3.7




######## UNNECESSARY FUNCTIONS ##################


# DEPRECATED
def extract_HyProtIDs(gff):
    HyProts = {}
    for entry in gff:
        try:
            attr = entry[8].split(';')
            ID = attr[0]                # here: keep "ID=" to make comparable to gff_content
            content = []
            for i in range(1,len(attr)):
                content.append(attr[i])
            prots.update({ID:content})
        except IndexError:
            pass
    return prots

#DEPRECATED
def get_HyProt_names(prot_dict):
    global HyProts_names
    regex = re.compile('hypothetical protein')  #search for HyProts
    if isinstance(prot_dict, dict):             # ensures loaded prot file is dictionary
        for prot in prot_dict.keys():
            # print(type(prot))
            # print(prot_dict[prot])
            if regex.search(str(prot_dict[prot])):
                HyProts_names.append(prot[3:])
            #     HyProts.append(prot)
    else:
        print(type(prot_dict))
        print("Proteins must be arranged in dictionary {ID:content} to extract hypothetical protein names.")
    # print(HyProts_names)

# DEPRECATED
def get_HyProt_seqs(input="/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/HyPro/test/prokka/chlamydia.ffn"):
    global HyProts, HyProts_names
    annotated = {}              # all annotated subsequences, dictionary
    header = ''
    seq = ''
    with open(input, 'r') as ffn:               # get all annotated header + seqs
        for line in ffn:
            if line.startswith('>') and re.search('hypothetical protein', line):    # header!
                elem = line.split(' ')
                if header != '' and seq != '':      # update only meaningful seqs
                    HyProts.update({header:seq})        # last sequence is complete, update the prots dictionary
                header = elem[0][1:]
                seq = ''
            elif re.match("[ACGT]",line):    # nt seq of hypothetical protein
                seq = seq + str(line)
        if header in HyProts_names: # if the last fasta seq is hypothetical protein, too
            elem = line.split(' ')
            print(header)
            HyProts.update({header:seq})
    print(len(HyProts.keys()))
    print(len(HyProts_names))
    # print(HyProts)
    assert len(HyProts) == len(HyProts_names)

    # for name in HyProts:
    #     try

                                                                             # exits program, if no sequence to write
    # os.system("grep -c '>'" + path + "/query.fasta")
    # assert int(os.system('grep -c '>'' + path + '/query.fasta')) == len(HyProts.keys())

# DEPRECATED
def is_gff(infile=in_gff):
    ''' Checks wether the file exists and extension is .gff'''
    EXT = get_input_info(infile)[2]
    if os.path.isfile(infile) and EXT == 'gff':
        print("Specified input is a file with 'gff' extension. Trying to load content...")
    else:
        print("Invalid input file. Please enter '--help' flag for information on tool usage.")
        exit()

# DEPRECATED
def is_outdir(out=out_dir):
    if os.path.isdir(out):
        print(f"Specified output directory is valid.")
    else:
        print("Invalid output directory. Please enter '--help' flag for information on script usage.")
        exit()

# DEPRECATED
#def format_notes(string, delimiter):
        # while len(elem) != 1:                           # asked for every new line
    #     if lines == 0:                          # add the qualifier after the first leading spaces
    #         l_string += '/notes=\"'
    #     else:
    #         l_string += (lspaces * ' ')
    #     lines += 1
    #     added = 0                           # for removal of elements after for loop ends
    #     for i in elem:
    #         nl_count = l_string.count('\n')                                           # newline counts; are counted in len() function, but not present in gbk format
    #         # print(f"newline counts \t {nl_count}")

    #         if (len(elem) - 1) != 0:
    #             length = len(l_string) + len(i) + len_del - nl_count      # + additional length of delimiter
    #         else:
    #             length = len(l_string) + len(i) + len_del + 1 - nl_count       # + additional quotes at the end of the string

    #         if length <= lines*80:                          # add the element if it does not exceed the 80 character mark
    #             added += 1
    #             info = i
    #             # elem.remove(i)
    #             l_string += info + ';'
    #         else:
    #             l_string +='\n'
    #             break                           # if the line would get longer than allowed break iteration over elements
    #     for i in range(0,added-1):                # erase added info from input string array
    #         elem.pop(i)
    # l_string = l_string.rstrip(';')
    # l_string += '\"\n'
    # print(l_string)
    # # exit()
    # return l_string
