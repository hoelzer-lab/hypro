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
hyprot_loc = {}         # recognizes the line of each hyprot via prokka feature ID
hyprot_content = {}            # recognizes prokka features of 'hypothetical proteins' via prokka feature ID
valid_db = ['uniprotkb']        # to be extended for novel db-types
id_scinames = {}               # feature IDs and target scientific name
################## ARGPARSE #####################
parser = argparse.ArgumentParser()
# Mode
parser.add_argument('-m', '--modus', dest='modus', metavar='[restricted, full]', choices=['restricted','full'], default = 'full', required = False, help="Modus of prokkaX to decide either full hyprot annotation or restricted (only full blanks). Valid arguments: 'full' and 'restricted'")
#Input-GFF
parser.add_argument('-i', '--input', dest='input', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the gff file, that shall be extended.")
# Output_dir
parser.add_argument('-o', '--output', dest='output', action='store', metavar='PATH', nargs=1, required=True, default='..', help="Specify PATH to a directory. prokkaX will generate all output to this.")
# DB-Type
parser.add_argument('-d', '--database', dest='db', action='store', metavar='STR', nargs=1, required=False, default='uniprotkb', help="Specifiy the target db to search in for annotation extension. Available options: 'uniprotkb'")
# location of mmseq2.sh:
parser.add_argument('-t', '--mmseqs2', dest='mmseq', action='store', metavar='PATH', nargs=1, required=True, help="Specify the path to the mmseqs2.sh. Obligatory for extension.")
# ALL or just BLANKS ?
#parser.add_argument('-m', '--modus', dest='modus', action='store', metavar='STR', nargs=1, default='blanks', help="Search all hypothetical proteins or just the blanks. default: blanks")



args = parser.parse_args()

if isinstance(args.modus, list):
    mode = args.modus[0]
elif isinstance(args.modus, str):
    mode = args.modus

print(f'Start prokkaX in {mode} mode')
in_gff = args.input[0]
out_dir = args.output[0]
if isinstance(args.db, list):
    db = args.db[0]
elif isinstance(args.db, str):
    db = args.db

# print(args.db)
# print(db)
# print(type(args.db))

ms = args.mmseq[0]


DIR = ''            # input directory, all prokka files should reside at
BN = ''             # basename of prokka output files
EXT = ''            # extension found in input
count = 0

################# FUNCTIONS ########################
def update_faa(in_gff, output):
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
    # print(faa_list)

    #write the updated faa file
    with open(output + '/output/' + BN + "_extended.faa", 'w') as faa:
        for entry in faa_list:
            faa.write(entry)

def update_ffn(in_gff, output):
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

                
    


# def update_fsa():
# def update_gbk():
# #and more:
# def update_sqn():
# def update_tbl():
# def update_tsv():
# def update_txt():

############ update_gff
def load_gff(input=in_gff):
    # is_outdir()
    # is_gff()
    row = 0
    hyprot_count = 0
    with open(input, 'r') as gff:
        for line in gff:
            line = line.rstrip('\n')
            attr = line.split('\t')
            gff_content.update({row:attr})      # ggf content with line ID as identifier
            if is_hyprot(attr, 'hypothetical protein'):
                hyprot_count += 1
                save_hyprot(attr, row)
            else:
                pass
            row += 1             # next line
    print("loading gff file finished.")
    print(f"Total Features:\t{row} ")
    print(f"Hyprot count:\t{hyprot_count}")
    # print(f"Hyprot_content:\t{len(hyprot_content.keys())}")
    # print(gff_content[len(gff_content)-1])
    print(f"Saved hyprots:\t{len(hyprot_content.keys())}")
    return gff_content, hyprot_content   # DEPRECATED

def is_hyprot(attr, regex = ''):
    '''Checks the attribute column for the regex "hypothetical protein" '''
    regex = re.compile(regex)
    try:
        if regex.search(attr[8]):
            return True
        else:
            return False
    except IndexError:
        pass

def save_hyprot(attr, row):
    '''get hyprot IDs, their attributes, and where in gff_content they can be found'''
    global hyprot_loc, hyprot_content, count, mode
    fields = attr[8].split(';')
    ID = fields[0][3:]
    content = []
    # print(f'{ID} {len(fields)} {fields[len(fields)-1]}')
    if mode == 'restricted':
        if len(fields) == 4:            # hypothetical proteins with no information at all
            for i in range(1, len(fields)):
                content.append(fields[i])
            # count += 1
            hyprot_content.update({ID:content})
            hyprot_loc.update({ID:row})
        else:
            pass 
        assert  len(hyprot_content.keys()) == len(hyprot_loc.keys())
    elif mode == 'full':
        if len(fields) == 4 or len(fields) == 5:            # hyprots with no info and partial info
            for i in range(1, len(fields)):
                content.append(fields[i])
            # count += 1
            hyprot_content.update({ID:content})
            hyprot_loc.update({ID:row})
        else:
            pass 
        assert  len(hyprot_content.keys()) == len(hyprot_loc.keys())
    # print(f"save_hyprots:\t{count}")
    return hyprot_content, hyprot_loc   # DEPRECATED

def query_fasta(infile = in_gff, output = out_dir):
    global hyprot_content, DIR, BN
    print("Build output directory 'prokkaX' in specified output location, if not existing.")
    print("Try to build query fasta...")
    ffn_file = DIR + '/' + BN + '.ffn'      # feature sequences
    if bool(hyprot_content):
        fasta = load_fasta(ffn_file)
        with open(output + "/query.fasta", 'w') as query:
            for ID in fasta.keys():
                # print(ID)
                if ID in hyprot_content.keys():
                    # print("hyprot")
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
    # print(len(hyprot_content.keys()))
    assert num == len(hyprot_content.keys()) 

def load_fasta(file):
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
            # if header in hyprots_names: # if the last fasta seq is hypothetical protein, too
            #     elem = line.split(' ')
            #     print(header)

# download a particular db
def download_db(output = out_dir, dbtype = db): #--------------------------TO BE EXTENDED -----------------------------------
    global valid_db
    weblink = ''
    if dbtype == valid_db[0]:    # if clause to decide which db to download 
        path = output.rstrip('/') + "/db/" + dbtype
        db_fasta = path + "/uniprot_sprot.fasta"
        db_target = path + "/target_db"
        if loaded(db_fasta, dbtype):
            print("uniprotkb fasta exists already. Skip download.")
        else:
            print("No db found. Download and index...")
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
    else:
        print("Error. Exit from download_db function...")
        exit()
    return path, db_fasta, db_target    #-----------------------------------------------------------------------------------------

def mmseq(dbfasta, dbtarget, output=out_dir, mmseq = ms):
    global hyprot_content, db, id_scinames 
    output = output.rstrip('/') + '/mmseq_output'
    # # execute mmseq2
    # os.system(f"{mmseq}  {output}  {dbfasta}  {dbtarget} {db}")
    
    # # fraction identified
    # hit_nums = str(subprocess.check_output("cut -f1 " + output + "/mmseq2_out_unique.tsv" + "| wc -l", shell=True))
    # print(f"Found hits for {int(hit_nums[2:-3])-1} / {len(hyprot_content.keys())} hypothetical proteins")

    # # create ID:annotation dictionary
    mmseqs_out = pd.read_csv(output + '/mmseq2_out_unique.tsv', sep='\t')

    #create dict from pandas df
    out_dict = mmseqs_out.to_dict('index')  # dict of dicts; {index:{row}}; row = {col-X:row-entry}
    # print(prots_dict.keys())

    # lookup sci_names and symbols of target IDs
    dict_scinames = get_names(output + '/mmseq2_out_unique.tsv')
    # print(dict_scinames)

    # save findings in hyprot_content (target info + scinames/symbols of target IDs)
    # print(hyprot_content)
    for index in out_dict.keys():           
        # print(out_dict[index]['query'])
        hyprot = out_dict[index]['query']
        db_ID = out_dict[index]['target']
        id_scinames.update({hyprot:dict_scinames[db_ID][1]})
        if hyprot in hyprot_content.keys():
            info = ''
            for attr in out_dict[index].keys():
                info = f"{info}aln_info_{attr}={out_dict[index][attr]};"
            info = f"{info}aln_info_symbol={dict_scinames[db_ID][0]};"      # add target symbol to info
            info = f"{info}aln_info_sciname={dict_scinames[db_ID][1]};"     # add target sci name to info
                # print(info)
                # print(out_dict[index][attr])
            # print(info)
            # print(hyprot_content[hyprot]) 
            hyprot_content[hyprot][0] = f"inference=mmseqs2 {db}"
            hyprot_content[hyprot][2] = f"product={out_dict[index]['target']}"
            hyprot_content[hyprot].append(info.rstrip(';'))
            # print(hyprot_content[hyprot]) 
            # print(hyprot_content[out_dict[index]['query']])
            # print(out_dict[index]['query'])
      
def update_gff(output=out_dir, delimiter = '\t'):
    global hyprot_content, hyprot_loc, gff_content, BN   
    output = output.rstrip('/') + "/output"
    for hyprot in hyprot_content.keys():        
        data = 'ID='
        data = f'{data}{hyprot}'
        # print(hyprot)
        # print(hyprot_content[hyprot])
        for elem in hyprot_content[hyprot]:
            data = f"{data};{elem}"
            # print(data)
        # data = data + '\''
        # print(data)
        loc = hyprot_loc[hyprot]
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
    global DIR, BN, EXT
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
    if os.path.isdir(output):
        print(f"Specified output directory is valid.")
    else:
        print("Invalid output directory. Please enter '--help' for information on script usage.")
        exit()
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
    if os.path.isfile(mmseq) and mmseq_name == "mmseqs2.sh":
        print("Found mmseqs2.sh.")
    else:
        print("mmseqs2.sh path specification is corrupted. Check the specified path and the file you referenced.")
        exit()
    
def loaded(dbfasta, db):
    global valid_db
    if db == valid_db[0]:
        if os.path.isfile(dbfasta):
            return True
        else:
            return False
    else:
        print("Unknown db type to be checked for loading status. (Message from 'loaded' method)")

def create_outdir(out_dir, dbtype):
    '''create the output path structure'''
    os.system(f'mkdir -p {out_dir}')
    os.system(f'mkdir -p {out_dir}/db')
    os.system(f'mkdir -p {out_dir}/db/{dbtype}')
    os.system(f'mkdir -p {out_dir}/mmseq_output/tmp')
    os.system(f'mkdir -p {out_dir}/output')
   
def get_names(table):   # table should be an input 
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
    ret = gene.querymany(ids, scopes='uniprot', fields='name,symbol', verbose= False)           # set verbose 'True' to get info about duplicates/missing values of name parsing
    sciname_dict = collect_scinames(ret)
    # print(sciname_dict)
    return sciname_dict
    # ret['query','name','symbol']
    # ret.to_csv('/mnt/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokkaX/mmseq_output/mmseq2_out_unique_names.tsv')
    # print(ret)
def collect_scinames(name_list):
    sciname_dict = {}               # {ID:[symbol, sciname]}
    for entry in name_list:
        try:
            sciname_dict.update({entry['query']:[entry['symbol'],entry['name']]})
        except: 
            sciname_dict.update({entry['query']:['hyprot','hypothetical protein']})
        # print(entry['query'])
    return sciname_dict


############# MAIN ######################

check_args()
create_outdir(out_dir, db)
load_gff()
query_fasta()
db_dir, dbfasta, dbtarget = download_db()
mmseq(dbfasta, dbtarget)
update_gff()
update_faa(in_gff, out_dir)
update_ffn(in_gff, out_dir)
# get_names('/mnt/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokkaX/mmseq_output/mmseq2_out_unique.tsv')


# TO DO:
    # object oriented - > atrributes: ID eC_number Name db_xref gene inference locus_tag product
    # design for different db
    # mmseq() - parametric input of target db info


# Requirements

# - prokka 1.14.0
# - mmseqs2 10.6d92c
# mygene module





######## UNNECESSARY FUNCTIONS ##################


# DEPRECATED
def extract_hyprotIDs(gff):     
    hyprots = {}
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
def get_hyprot_names(prot_dict):
    global hyprots_names
    regex = re.compile('hypothetical protein')  #search for hyprots
    if isinstance(prot_dict, dict):             # ensures loaded prot file is dictionary
        for prot in prot_dict.keys():
            # print(type(prot))
            # print(prot_dict[prot])
            if regex.search(str(prot_dict[prot])):
                hyprots_names.append(prot[3:])   
            #     hyprots.append(prot)
    else:
        print(type(prot_dict))
        print("Proteins must be arranged in dictionary {ID:content} to extract hypothetical protein names.")
    # print(hyprots_names)

# DEPRECATED
def get_hyprot_seqs(input="/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokka/chlamydia.ffn"):
    global hyprots, hyprots_names
    annotated = {}              # all annotated subsequences, dictionary
    header = ''
    seq = ''
    with open(input, 'r') as ffn:               # get all annotated header + seqs
        for line in ffn:
            if line.startswith('>') and re.search('hypothetical protein', line):    # header!
                elem = line.split(' ')
                if header != '' and seq != '':      # update only meaningful seqs
                    hyprots.update({header:seq})        # last sequence is complete, update the prots dictionary
                header = elem[0][1:]
                seq = '' 
            elif re.match("[ACGT]",line):    # nt seq of hypothetical protein
                seq = seq + str(line)
        if header in hyprots_names: # if the last fasta seq is hypothetical protein, too
            elem = line.split(' ')
            print(header)
            hyprots.update({header:seq})
    print(len(hyprots.keys()))
    print(len(hyprots_names))
    # print(hyprots)
    assert len(hyprots) == len(hyprots_names)

    # for name in hyprots:
    #     try

                                                                             # exits program, if no sequence to write
    # os.system("grep -c '>'" + path + "/query.fasta")
    # assert int(os.system('grep -c '>'' + path + '/query.fasta')) == len(hyprots.keys())

# DEPRECATED
def is_gff(infile=in_gff):
    ''' Checks wether the file exists and extension is .gff'''
    EXT = get_input_info(infile)[2]
    if os.path.isfile(infile) and EXT == 'gff':
        print("Specified input is a file with 'gff' extension. Trying to load content...")
    else:
        print("Invalid input file. Please enter '--help' flag for information on script usage.")
        exit()

# DEPRECATED
def is_outdir(out=out_dir):
    if os.path.isdir(out):
        print(f"Specified output directory is valid.")
    else:
        print("Invalid output directory. Please enter '--help' flag for information on script usage.")
        exit()