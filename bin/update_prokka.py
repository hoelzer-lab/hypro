#!/usr/bin/env python3

####### PACKAGES ##############################################################
import re, os, sys, subprocess, argparse
import mygene
import pandas as pd
import json

###### ############GLOBAL VAR ################################################
valid_db = ['uniprotkb','uniref100', 'uniref90', 'uniref50', 'pdb']        # to be extended for novel db-types
id_scinames = {}               # feature IDs and target scientific name

################## ARGPARSE ####################################################
parser = argparse.ArgumentParser()
# Mode
parser.add_argument('-m', '--modus', dest='modus', metavar=['restricted', 'full'], choices=['restricted','full'], default = 'full', required = False, help="Modus of HyPro to decide either for an all hypothetical protein annotation or restricted (only full blanks with no partial annotation). Valid arguments: 'full' and 'restricted'")
#Input-FFN
parser.add_argument('-i_ffn', '--input_ffn', dest='input_ffn', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the ffn file, that shall be extended.")
#Input-FAA
parser.add_argument('-i_faa', '--input_faa', dest='input_faa', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the faa file, that shall be extended.")
#Input-FAA
parser.add_argument('-i_gbk', '--input_gbk', dest='input_gbk', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the gbk file, that shall be extended.")
#Input from mmseqs2 Output
parser.add_argument('-ms', '--mmseqs2', dest='mmseqs2', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the *unique.tsv file from mmseqs2 output")
#Hyprot_loc
parser.add_argument('-hl', '--hyprot_loc', dest='hyprot_loc', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the hyprot_loc dict")
#Hyprot_content
parser.add_argument('-hc', '--hyprot_content', dest='hyprot_content', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the hyprot_content dict")
#gff_content
parser.add_argument('-gffc', '--gff_content', dest='gff_content', action='store', metavar='PATH', nargs=1, required=True, help="Specify PATH to the gff_content dict")
# Output_dir
parser.add_argument('-o', '--output', dest='output', action='store', metavar='PATH', nargs=1, required=True, default='..', help="Specify PATH to a directory. HyPro will generate all output to this.")
# Prokka files basename
parser.add_argument('-n', '--name', dest='name', action='store', metavar='PATH', nargs=1, required=True, help="Specify basename of prokka output files that will be extended.")
# DB-Type
parser.add_argument('-d', '--database', dest='db', action='store', metavar='STR', nargs=1, required=False, default='uniprotkb', help="Specify the target db to search in for annotation extension. Available options: 'uniprotkb', 'uniref100', 'uniref90', 'uniref50', 'pdb' [uniprotkb]")
args = parser.parse_args()

if isinstance(args.mmseqs2, list):
    ms = args.mmseqs2[0]
elif isinstance(args.mmseqs2, str):
    ms = args.mmseqs2

if isinstance(args.input_ffn, list):
    in_ffn = args.input_ffn[0]
elif isinstance(args.input_ffn, str):
    in_ffn = args.input_ffn

if isinstance(args.input_faa, list):
    in_faa = args.input_faa[0]
elif isinstance(args.input_faa, str):
    in_faa = args.input_faa

if isinstance(args.input_gbk, list):
    in_gbk = args.input_gbk[0]
elif isinstance(args.input_gbk, str):
    in_gbk = args.input_gbk

if isinstance(args.modus, list):
    mode = args.modus[0]
elif isinstance(args.modus, str):
    mode = args.modus

if isinstance(args.name, list):
    name = args.name[0]
elif isinstance(args.name, str):
    name = args.name

if isinstance(args.gff_content, list):
    gffc = args.gff_content[0]
elif isinstance(args.gff_content, str):
    gffc = args.gff_content

if isinstance(args.hyprot_loc, list):
    hl = args.hyprot_loc[0]
elif isinstance(args.hyprot_loc, str):
    hl = args.hyprot_loc

if isinstance(args.hyprot_content, list):
    hc = args.hyprot_content[0]
elif isinstance(args.hyprot_content, str):
    hc = args.hyprot_content

with open(hl, 'r') as f:
    HyProt_loc = json.load(f)         # recognizes the line of each HyProt via prokka feature ID
with open(hc, 'r') as f:
    HyProt_content = json.load(f)            # recognizes prokka features of 'hypothetical proteins' via prokka feature ID
with open(gffc, 'r') as f:
    gff_content = json.load(f)        # gff content per row number

if isinstance(args.db, list):
    db = args.db[0]
elif isinstance(args.db, str):
    db = args.db

out_dir = args.output[0]

################# FUNCTIONS ##################################################
def update_faa(output=out_dir, infile=in_faa): #former 2nd arg: in_gff
    global DIR, BN, id_scinames
    faa_list = []
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
    with open(output + '/' + name + "_extended.faa", 'w') as faa:
        print(f"Write extended faa:\t{name}_extended.faa")
        for entry in faa_list:
            faa.write(entry)

def update_ffn(output=out_dir, infile=in_ffn): #former 2nd arg: in_gff
    global DIR, BN, id_scinames
    ffn_list = []
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
    with open(output + '/' + name + "_extended.ffn", 'w') as ffn:
        print(f"Write extended ffn:\t{name}_extended.ffn")
        for entry in ffn_list:
            ffn.write(entry)

############ extend gbk
def extend_gbk(output, id_alninfo, infile=in_gbk):
    '''extend the information in prokka output gbk with mmseq2 homology findings - based on the description of genbank file structure of NCBI'''
    global id_scinames
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
    with open(output + "/" + name + "_extended.gbk", 'w') as gbk_out:
        print(f"Write extended gbk:\t{name}_extended.gbk")
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
def get_mmseq_output(mmseqs2_out=ms):
    '''
    Search hypothetical proteins from prokka output in a selected database using mmseqs2.
    Extract, format and write new annotation to output id_infos, update previous HyProt_content.
    '''
    global HyProt_content, db, id_scinames
    id_infos = {}                                           # additional information found  with mmseq; saved with ID as key
    # fraction identified
    hit_nums = str(subprocess.check_output(f"cut -f1 {ms} | wc -l", shell=True))
    print(f"The homology search assigns a function to {int(hit_nums[2:-3])-1} / {len(HyProt_content.keys())} hypothetical proteins\n\t(this includes also hits on hypothetical proteins in the database)")

    # load annotation pandas dataframe
    mmseqs_out = pd.read_csv(mmseqs2_out , sep='\t')

    #create dict from pandas df - ID:annotation dictionary
    out_dict = mmseqs_out.to_dict('index')  # dict of dicts; {index:{row}}; row = {col-X:row-entry}
    # if out_dict is unref50/90/100:
    if db == valid_db[1] or db == valid_db[2] or db == valid_db[3]:
        for num in out_dict.keys():
            # out_dict[num]["target"]
            value = out_dict[num]
            target_split = out_dict[num]['target'].split('_')
            out_dict[num].update({'target':target_split[1]})  # ID form "Uniref[50,90,100]_uniprotID" replaced by uniprotID


    # lookup sci_names and symbols of target IDs
    dict_scinames = get_names(mmseqs2_out, db)

    #print(f'SciNames of TargetIDs:\n{dict_scinames}')
    # save findings in HyProt_content (target info + scinames/symbols of target IDs)
    for index in out_dict.keys():
        HyProt = out_dict[index]['query']
        db_ID = out_dict[index]['target']
        id_scinames.update({HyProt:dict_scinames[db_ID][1]})                # save scinames under the HyProt prokka IDs
        if HyProt in HyProt_content.keys():
            info = ''
            for attr in out_dict[index].keys():
                info = f"{info}HyPro_{attr}={out_dict[index][attr]};"
            info = f"{info}HyPro_symbol={dict_scinames[db_ID][0]};"      # add target symbol to info
            info = f"{info}HyPro_sciname={dict_scinames[db_ID][1]};"     # add target sci name to info
            id_infos[HyProt] = [info]
            id_infos[HyProt].append(out_dict[index]['target'])
            HyProt_content[HyProt][0] = f"inference=mmseqs2 {db}"
            HyProt_content[HyProt][2] = f"product={out_dict[index]['target']}"
            HyProt_content[HyProt].append(info.rstrip(';'))
    counts = count_HyProts(id_scinames)
    print(f"Of those: {counts} / {int(hit_nums[2:-3])-1} are hypothetical proteins.")
    return id_infos

def update_gff(output=out_dir, delimiter = '\t'):
    global HyProt_content, HyProt_loc, gff_content
    for HyProt in HyProt_content.keys():
        data = 'ID='
        data = f'{data}{HyProt}'
        for elem in HyProt_content[HyProt]:
            data = f"{data};{elem}"

        loc = HyProt_loc[HyProt]
        gff_content[str(loc)][8] = data
    # write to gff
    with open(output +'/' + name + "_extended.gff", 'w') as gff_up:
        print(f"Write extended gff:\t{name}_extended.gff")
        for key in gff_content:
            line = ''
            for elem in gff_content[key]:
                line = line + str(elem) + '\t'
            # break
            gff_up.write(str(line.rstrip('\t')) + '\n')

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



################################  MAIN  ########################################
id_alninfo = get_mmseq_output()
update_gff()
update_faa()
update_ffn()
extend_gbk(out_dir, id_alninfo)
