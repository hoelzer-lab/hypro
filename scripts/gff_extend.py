#!/usr/bin/env python3

import re, os, sys, subprocess
# import mygene
# from biothings_client import get_client
import pandas as pd
from Bio import SeqIO

# import wget
from collections import OrderedDict
db = 'uniprotkb'          # shall be paramtric to choose between different databases...

gff_content = {}        # gff content per row number
hyprot_loc = {}         # recognizes the line of each hyprot via prokka feature ID
hyprot_content = {}            # recognizes prokka features of 'hypothetical proteins' via prokka feature ID

def load_gff(input="/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokka/chlamydia_test.gff"):
    row = 0
    with open(input, 'r') as gff:
        for line in gff:
            line = line.rstrip('\n')
            attr = line.split('\t')
            gff_content.update({row:attr})      # ggf content with line ID as identifier
            if is_hyprot(attr, 'hypothetical protein'):
                save_hyprot(attr, row)
            else:
                pass
            row += 1             # next line
    # print(gff_content[len(gff_content)-1])
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
    global hyprot_loc, hyprot_content
    fields = attr[8].split(';')
    ID = fields[0][3:]
    content = []
    if len(fields) == 4:            # only hypothetical proteins with no information at all! 
        for i in range(1, len(fields)):
            content.append(fields[i])
        hyprot_content.update({ID:content})
        hyprot_loc.update({ID:row})
    else:
        pass 
    assert  len(hyprot_content.keys()) == len(hyprot_loc.keys())
    return hyprot_content, hyprot_loc   # DEPRECATED

def query_fasta(input= "/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokka/", output='/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x/query.fasta'):
    global hyprot_content
    if bool(hyprot_content):
        fasta = load_fasta(input + '/chlamydia.ffn')
        with open(output, 'w') as query:
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
        print("No query seq to write to fasta. Exiting extension pipeline...")      
        exit()   
    num =  int(subprocess.check_output("grep -c '>' " + output, shell=True))
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
    print(len(fasta.keys()))
    return fasta
    # print(fasta)        
            # if header in hyprots_names: # if the last fasta seq is hypothetical protein, too
            #     elem = line.split(' ')
            #     print(header)



# download a particular db
def download_db(dir = "/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x"):
    path = dir + '/db/uniprotkb'
    # print(path)
    os.system('mkdir -p' + " " + path) # uniprotkb - make parametric
    os.chdir(path)
    os.system("wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz")    
    os.system('gunzip ' + path + "/uniprot_sprot.fasta.gz")


def mmseq(path='/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x'):
    global hyprot_content, db
    # execute mmseq2
    os.system('/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/scripts/mmseq2.sh /data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x')
    
    # fraction identified
    hit_nums = str(subprocess.check_output('cut -f1 /data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x/mmseq2_out_unique.tsv | wc -l',shell=True))
    print(f"Found hits for {int(hit_nums[2:-3])-1} / {len(hyprot_content.keys())} hypothetical proteins")
    
    
    # create ID:annotation dictionary
    mmseqs_out = pd.read_csv('/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x/mmseq2_out_unique.tsv', sep='\t')

    out_dict = mmseqs_out.to_dict('index')  # dict of dicts; {index:{row}}; row = {col-X:row-entry}
    # print(prots_dict.keys())
    
    for index in out_dict.keys():           # aln
        # print(out_dict[index]['query'])
        hyprot = out_dict[index]['query']
        if hyprot in hyprot_content.keys():
            info = ''
            for attr in out_dict[index].keys():
                info = f"{info}aln_info_{attr}={out_dict[index][attr]};"
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
            
def update_gff(output="/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/prokka/chlamydia_test_extend.gff", delimiter = '\t'):
    global hyprot_content, hyprot_loc, gff_content        
    for hyprot in hyprot_content.keys():        # keine
        data = 'ID='
        data = f'{data}{hyprot}'
        # print(hyprot)
        # print(hyprot_content[hyprot])
        for elem in hyprot_content[hyprot]:
            data = f'{data};{elem}'
            # print(data)
        # data = data + '\''
        # print(data)
        loc = hyprot_loc[hyprot]
        # print(gff_content[loc])
        gff_content[loc][8] = data
        # print(gff_content[loc])
    # write to gff
    with open(output, 'w') as gff_up:
        for key in gff_content:
            line = ''
            # print(gff_content[key])
            for elem in gff_content[key]:
                line = line + str(elem) + '\t'
                # print(line)
            # break
            #     print(line)
            gff_up.write(str(line.rstrip('\t')) + '\n')


            
   
    



# probleme: ambigious symbols, 1/3 IDs not found
def get_names(table):   # table should be an input 
    gene_ids = []
    translate = {}
    ids = []
    with open(table,'r') as csv:
        for line in csv:
            elem = line.split('\t')
            gene = mygene.MyGeneInfo()
            ids.append(elem[1])
            # print(gene.getgene('uniport:P24941','name,symbol'))
    ret = gene.querymany(ids, scopes='uniprot', fields='name,symbol', as_dataframe=True)
    ret['query','name','symbol']
    # ret.to_csv('/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x/notfound.csv')
    # print(ret)




load_gff()
query_fasta()
download_db()
mmseq()
update_gff()
# get_names('/data/mahlzeitlocal/projects/ma_neander_assembly/hiwi/prokkaX/test/x/mmseq2_out_unique.tsv')




# TO DO:
    # object oriented - > atrributes: ID eC_number Name db_xref gene inference locus_tag product
    # argparse
    # design for different db
    # ALL or only the hard cases
    # 


# Requirements

# - prokka 1.14.0
# - mmseqs2 10.6d92c




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

