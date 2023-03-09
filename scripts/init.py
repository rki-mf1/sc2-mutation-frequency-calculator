#!/usr/bin/env python3

# Research Internship FU Berlin
# RKI MF1
# Ashkan Ghassemi
# Variant specific PCR finder

from sys import exit
from os import path
import pandas as pd
from Bio.Seq import MutableSeq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter
import itertools
import numpy as np
import pandas as pd 
import yaml 
import os.path 
import argparse
import re
from pathlib import Path
from functools import reduce

########## input ############

###Helperfunction###
def txt_to_string(mutation_file):  # 'lineages/'aa', 'lineages.txt'/'aa_changes.txt'
    ''' converts a txt file (usually as args input) in a list
    input: either lineages or aa_mutations txt file
    output: list of either lineages without parent information or single aa mutations
    '''
    if "lineages" in mutation_file:
        # 2) lineages.txt
        # have to consider VOC/VOI/...
        lineage_txt = "src/lineages.txt"
        list_lineages_txt = Path(
            lineage_txt).read_text().replace('\n', '').split(",")
        pangolin_lineages = [i.split(",")[1] for i in Path(
            lineage_txt).read_text().split("\n")[:-1]]
        mutation_string = ' '.join(pangolin_lineages)
    if "aa" in mutation_file:
        # read .txt in list of strings
        # 1) aa_changes.txt
        # have to consider single/combinations of mutations
        aa_change_txt = "src/aa_changes.txt"
        list_aa_txt = Path(aa_change_txt).read_text().replace(
            '\n', ',').split(",")[:-1]
        list_comb_aa = [comb for comb in list_aa_txt if " " in comb]
        list_single_aa = [
            single for single in list_aa_txt if single not in list_comb_aa]
        list_single_aa = ' '.join(list_single_aa)
        mutation_string = [list_single_aa, list_comb_aa]
    return mutation_string

#for cut_off to be between 0.0 and 1.0
def restricted_float(x):
    ''' restricts input value of cut-off to be between 0.0 and 1.0
    input: cut-off value (args) 
    output: error if value is out of range or the none
    '''
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError(
            "%r not a floating-point literal" % (x,))

    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x

########## main functionality ############

def init_num_lineages(sample_column, mutation_profile):
    ''' init to calculate absolut occurences of lineages  
    input: 
        sample column: either lineages (default) or accession (covsonar check)
        mutation profile: dataframe from covsonar output  
    output: 
        dna/aa profile: dataframe from covsonar output (relevant columns)
        num lineage: dict with lineages and their respective occurrences in mutation profile
    '''
    df_mutation_profile = pd.read_table(mutation_profile)
    df_dna_aa_profile = df_mutation_profile[[
        sample_column, 'dna_profile', 'aa_profile']]
    # compute absolute values of mutations
    df_num_lineage = df_mutation_profile[sample_column].value_counts(
    ).rename_axis('lineage').reset_index(name='counts')
    dict_num_lineage = pd.Series(
        df_num_lineage.counts.values, index=df_num_lineage.lineage).to_dict()
    return df_dna_aa_profile, dict_num_lineage

def lineage_mutation_frequency(mutation_type, df_mutation_profile, lineages, num_lineages):
    ''' calculates relative mutation numbers based on init and options 
    input: 
        mutation_type: "dna_profile" or "aa_profile"
        df_mutation_profile: df consists of lineage, dna_profile, aa_profile
        lineages: list of lineages
        num_lineages: dict of lineages and their number of occurence
    output: dict with lineages and their respective occurrences in mutation profile
    '''
    list_df_lin_frequency = []
    for lineage in lineages:
        df_lin_frequency = pd.DataFrame()
        df_lineage_subprofile = df_mutation_profile.loc[df_mutation_profile['lineage'] == lineage]
        if not df_lineage_subprofile.empty:
            list_mutations = list(itertools.chain(
                *[i.split() for i in df_lineage_subprofile[mutation_type].unique().tolist()]))
            df_num_mutations = pd.Series(list_mutations).value_counts(
                ).rename_axis('mutation').reset_index(name='counts')
            dict_num_mutations = pd.Series(
                df_num_mutations.counts.values, index=df_num_mutations.mutation).to_dict()
            dict_num_lineage_mutations = {lineage: dict_num_mutations}
            df_lin_frequency[lineage] = pd.DataFrame.from_dict(
                dict_num_lineage_mutations[lineage], orient='index') / num_lineages[lineage]
            list_df_lin_frequency.append(df_lin_frequency)
    merged_df = reduce(lambda left, right: pd.merge(
        left, right, left_index=True, right_index=True, how='outer'), list_df_lin_frequency)
    return merged_df


########## optional parameters ############
## Matrix

def create_matrix(args):
    return

def create_heatmap_matrix():
    return

## Consensus

def snp_consensus(fasta, mutation):
    ''' definition of SNP mutation  
    input: 
        fasta: fasta sequence (Wuhan reference genome) 
        mutation: mutation in format ref_base:ref_position:alt_base
    output: altered fasta sequence
    '''
    consensus_seq_snp = fasta
    match = re.match(
        r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
    if match:
        items = list(match.groups())
        consensus_seq_snp = consensus_seq_snp[:int(items[1])-1] + \
            items[2] + consensus_seq_snp[int(items[1]):]
        fasta = consensus_seq_snp
    return fasta

def insert_consensus(fasta, mutation):
    ''' definition of Insertion mutation
    input:
        fasta: fasta sequence (altered fasta sequence after deletions)
        mutation: mutation in format ref_base:ref_position:(many)alt_bases
    output: altered fasta sequence
    '''
    consenus_seq_insert = fasta
    match = re.match(
        r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
    if match:
        items = list(match.groups())
        consenus_seq_insert = consenus_seq_insert[:int(items[1])-1] + \
            items[2] + consenus_seq_insert[int(items[1]):]
    fasta = consenus_seq_insert
    return fasta

#DELETION
def del_consensus(fasta, mutation):
    ''' definition of Deletion mutation
    input:
        fasta: fasta sequence (altered sequence after SNPs)
        mutation: mutation in format del:ref_position:number of bases
    output: altered fasta sequence
    '''
    consensus_seq_del = fasta
    deletion = mutation.split(":")
    consensus_seq_del = consensus_seq_del[:int(
        deletion[1])-1] + "?"*int(deletion[2]) + consensus_seq_del[int(deletion[1])+int(deletion[2])-1:]
    fasta = consensus_seq_del
    return fasta

def create_consensus(infile, lineage_dna_frequency, threshold):
    ''' creates many fasta for respective lineages using SNP, Deletions and Insertions 
    input: 
        infile: path for Wuhan ref seq
        lineage dna frequency: dataframe with frequency of mutations from their respective lineages
        threshold: cut-off for frequency (default: 0.75) 
    output: many fasta 
    '''  
    with open(infile) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            id, record = str(record.id), str(record.seq)
    consensus_id, ref_seq = MutableSeq(id), MutableSeq(record)

    df_consensus_mutations = pd.DataFrame(
        columns=['lineage', 'dna_mutations'])

    for column in lineage_dna_frequency.columns:
        list_mutations = lineage_dna_frequency.index[lineage_dna_frequency[column]
                                                            >= threshold].tolist()
        if list_mutations: #later important for checking consensus seqs
            df_consensus_mutations = df_consensus_mutations.append(
                {'lineage': column, 'dna_mutations': list_mutations}, ignore_index=True)
        
        dict_consensus_mutations = {column: list_mutations}
        for lineage, list_mutations in dict_consensus_mutations.items():
            snp, insertion, deletion = ([] for i in range(3))
            consensus_id, consensus_seq = lineage, ref_seq

            file_out = f'output/consensus/sc2_consensus_{lineage}.fasta'
            for mutation in list_mutations:
                # deletion
                if ":" in mutation:
                    deletion.append(mutation)
                else:
                    match = re.match(r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
                    if match:
                        items = list(match.groups())
                        # insertion
                        if len(items[2]) != 1:
                            insertion.append(mutation)
                        else:
                            # SNP
                            snp.append(mutation)
            if snp:
                for mutation in snp:  
                    print(mutation)
                    consensus_seq_snp = consensus_seq
                    match = re.match(
                        r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
                    if match:
                        items = list(match.groups())
                        ref_position = int(items[1])
                        ref_base = consensus_seq[int(items[1])]
                        ref_mutation_base = items[0]
                        print(f"Ref Seq at position {ref_position} is {ref_base}, ref base from mutation is {ref_mutation_base}")
                        # warning message
                        '''
                        if consensus_seq_snp[int(items[1])] != items[0]:
                            print(
                                "Warning: reference nucleotide doesnt match the snp mutation")
                        '''
                        consensus_seq_snp = consensus_seq_snp[:int(items[1])-1] + \
                            items[2] + consensus_seq_snp[int(items[1]):]
                        consensus_seq = consensus_seq_snp
            if deletion:
                for mutation in deletion:
                    consensus_seq_del = consensus_seq
                    deletion = mutation.split(":")
                    consensus_seq_del = consensus_seq_del[:int(
                        deletion[1])-1] + "?"*int(deletion[2]) + consensus_seq_del[int(deletion[1])+int(deletion[2])-1:]
                    consensus_seq = consensus_seq_del
            if insertion:
                for mutation in insertion:
                    consenus_seq_insert = consensus_seq
                    match = re.match(
                        r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
                    if match:
                        items = list(match.groups())
                        #warning message
                        if consenus_seq_insert[int(items[1])] != items[0]:
                            print(
                                "Warning: reference nucleotide doesnt match the insert mutation")
                        consenus_seq_insert = consenus_seq_insert[:int(items[1])-1] + \
                            items[2] + consenus_seq_insert[int(items[1]):]
                    consensus_seq = consenus_seq_insert

                        
            consensus_seq = MutableSeq(str(consensus_seq).replace("?", ""))
            seq_record = SeqRecord(
                    consensus_seq, id=consensus_id, description="NC_045512.2 Consensus Sequence")

            with open(file_out, "w") as f:
                SeqIO.write(seq_record, f, "fasta")
            print(f"Consensus file for {consensus_id} is created in {file_out}")
    print("Consensus sequences can now be analysed using MSA and phylogenetic trees")
    return df_consensus_mutations

def create_merged_fasta():
    ''' create merged fasta from many fasta 
    input: many fasta 
    output: merged fasta
    '''
    return

def consensus_mutations(lineages, df_mutation_profile):
    ''' creates dataframe of lineages and there mutations for consensus check 
    input: list of lineages 
    output: dataframe of lineages and there mutations (accession instead of lineages)
    '''
    df_consensus_mutations = pd.DataFrame(
        columns=['lineage', 'dna_mutations'])
    for lineage in lineages:
        df_lineage_subprofile = df_mutation_profile.loc[
            df_mutation_profile['accession'] == lineage]
        if not df_lineage_subprofile.isnull().values.any():
            list_mutations = list(itertools.chain(
                *[i.split() for i in df_lineage_subprofile['dna_profile'].unique().tolist()]))
            df_consensus_mutations = df_consensus_mutations.append(
                {'lineage': lineage, 'dna_mutations': list_mutations}, ignore_index=True)
    return df_consensus_mutations

def compare_lists(l1, l2):
    ''' compares two lists 
    input: two unsorted lists 
    output: boolean 
    '''
    return set(l1) == set(l2)



# entry point
def main():
    # Container to hold the arguments
    parser = argparse.ArgumentParser(description='VOC PCR Finder')
    # input
    # consensus related, --sequences ?
    parser.add_argument('-tsv', '--mutation_profile', metavar='', required=False,
                        help='choose mutation_profile (TSV) generated by covsonar from input directory')
    parser.add_argument('-db', '--database', metavar='', required=False,
                        help='database including virus variant whole genomes send to RKI daily')
    parser.add_argument('-date', '--date_range', metavar='', required=False,
                        help='optional: statistic will be computed for sequences of given time range (y-m-d:y-m-d), e.g.: 2021-05-01: 2021-05-07')
    parser.add_argument('-lin', '--lineages', metavar='', required=False,
                        help='optional: choose specific VOI/VOC lineages which are set in TXT file')
    parser.add_argument('-aa', '--aa_mutations', metavar='', required=False,
                        help='optional: choose specific aa mutations which are set in TXT file')
    parser.add_argument('-m', '--matrix', required=False, action='store_true',
                        help='optional: create matrix from VOI/VOC lineages and there mutation frequencies')
    parser.add_argument('-map', '--heatmap', required=False, action='store_true',
                        help='optional: create heatmap visualization from ceated matrix')
    parser.add_argument('-con', '--consensus', required=False, action='store_true',
                        help='optional: will create consensus sequences based on given lineages')
    parser.add_argument('-check', '--consensus_check', required=False, action='store_true',
                        help='optional: alterate path to check whether the created consensus seqs have the right mutations')
    parser.add_argument('-sig', '--signature', required=False, action='store_true',
                        help='optional: will give all signature mutation for given lineages')
    # https://bioinformaticshome.com/db/collection/phylogenetics
    parser.add_argument('-vis', '--visualization', metavar='', required=False,
                        help='optional: visualization of consensus sequences')
    parser.add_argument('-cut', '--cut_off', metavar='', required=False, type=restricted_float,
                        help='optional: set threshold for aa mutations to create matrix or consensus')
    args = parser.parse_args()

    #load yml file 
    with open("config/config.yml", "r") as ymlfile:
        config = yaml.full_load(ymlfile)
    
    # Step1: check if covsonar mutation profile (output tsv) is given
    #AA Mutation als Fall für tsv generation DATAFRAME MERGE ÄNDERN AUF OUTTER
    if not args.mutation_profile:  
        print("Mutation profile from Covsonar doesnt exist: \n"
              "Covsonar running ...")
        if args.lineages or args.aa_mutations:
            if args.lineages:
                lineages_string = txt_to_string(args.lineages)
                outfile = 'input/mutation_lineages_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match --lineage {lineages_string} --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}. \n"
                    "Choosen Lineages: \n"
                    f"{lineages_string}")
            if args.aa_mutations: #für jede mutation (single als auch comb) muss covsonar einzeln ausgeführt werden, freq berechnung für single/comb dann einzeln basierend auf merged tsv
                aa_single_changes, aa_comb_changes = txt_to_string(args.aa_mutations)[0].replace(
                    " ", " " + '-i' + " "), txt_to_string(args.aa_mutations)[1]
                outfile = 'input/mutation_aa_single_profile.tsv'
                #comb_outfiles = 'input/combinations/mutation_aa{counter}_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match -i {aa_single_changes} --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}. \n"
                    "Choosen Amino acid Changes: \n"
                    f"{aa_single_changes}")
        elif args.lineages and args.aa_mutations:
            lineages_string = txt_to_string(args.lineages)
            aa_single_changes = txt_to_string(args.aa_mutations)[0].replace(
                " ", " " + '-i' + " ")
            outfile = 'input/mutation_aa_lineages_profile.tsv'
            os.system(
                f"python3 covsonar/sonar.py match -i {aa_single_changes} --lineage {lineages_string} --date {args.date_range} --db {args.database} --tsv > {outfile}")
            print(f"Covsonar created mutation profile in {outfile}. \n"
                "Choosen Amino acid Changes: \n"
                f"{aa_single_changes} \n"
                "Choosen Lineages: \n"
                f"{lineages_string}")
        else:
            if not args.date_range or not args.database:
                print("A date range and database should be available for Covsonar to create mutation profile")
            else:
                outfile = 'input/mutation_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}.")
    
    else: # db and date not more relevant;
        if os.path.exists(args.mutation_profile):
            print("A mutation profile (covsonar output) exists")
            print("Read " + str(args.mutation_profile))

            df_dna_aa_profile, dict_num_lineage = init_num_lineages(
                'lineage', args.mutation_profile)

            print("Compute statistics ...")
            
            # compute absolute values of mutations (dict)
            if args.lineages:  # .txt has to be string
                lineage_list = txt_to_string(args.lineages).split(' ')
            else:
                lineage_list = list(dict_num_lineage.keys())
            

            if args.matrix: 
                print("Matrix will be created..")
                df_lineage_mutation_frequency = lineage_mutation_frequency(
                    "aa_profile", df_dna_aa_profile, lineage_list, dict_num_lineage)
                
                if not args.cut_off:  # rowwise
                    print("no cutoff")
                    if not args.aa_mutations:
                        print("and no aa mutations")
                        df_matrix = df_lineage_mutation_frequency
                        print(df_matrix)
                    else:
                        print("but aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                        print(df_matrix)
                else: #cutoff
                    print("cutoff is selected")
                    cut_off = args.cut_off
                    if not args.aa_mutations:
                        print("but no aa mutations")
                        df_matrix = df_lineage_mutation_frequency.fillna(
                            0)[(df_lineage_mutation_frequency.fillna(0) >= cut_off).any(axis=1)]
                        print(df_matrix)
                    else: #cutoff + aa 
                        print("and aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix_aa = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                        df_matrix = df_matrix_aa.fillna(0)[(df_matrix_aa.fillna(0) >= args.cut_off).any(axis=1)]
                        print(df_matrix)
                df_matrix.to_excel("output/matrix/mutation_frequency.xlsx")
            
            elif args.consensus:
                print("Mutation Frequency per lineage wil be computed..")

                df_lineage_mutation_frequency = lineage_mutation_frequency(
                    "dna_profile", df_dna_aa_profile, lineage_list, dict_num_lineage)
                
                print("creating consensus sequences...")
                file_in = 'data/NC_045512.2.fasta.txt'
                
                if not args.cut_off:  # cut-off columnwise
                    print("cut-off = 0.75 since no cut-off is selected")
                    cut_off = 0.75
                else:
                    cut_off = args.cut_off
                    print(f"cut-off is set to {cut_off}")

                df_created_consensus_mutations = create_consensus(file_in, df_lineage_mutation_frequency, cut_off)
                
                # Analyse Consensus Sequences
                fasta_dir = "output/consensus/"
                records = []

                for filename in os.listdir(fasta_dir):
                    if filename.endswith(".fasta") or filename.endswith(".fa"):
                        filepath = os.path.join(fasta_dir, filename)
                        for record in SeqIO.parse(filepath, "fasta"):
                            records.append(record)

                # Write the merged fasta sequences to a new file
                with open("output/consensus/merged.fasta", "w") as f:
                    SeqIO.write(records, f, "fasta")
                print(
                    "Merged consensus file is created in output/consensus/merged.fasta")
                
                if args.consensus_check:
                    print("consensus seq check if the mutations are right")
                    print("Created Consenus will be checked via Covsonar..")
                    print("Mutation profiles will be generated")

                    os.system('touch data/consensus.db')
                    os.system(
                        'python3 covsonar/sonar.py add -f output/consensus/merged.fasta --db data/consensus.db --cpus 8 --noprogress')
                    os.system('python3 covsonar/sonar.py match --db data/consensus.db --tsv > input/consensus_mutation_profile.tsv')
                    
                    print("Mutations will be checked...")
                    df_dna_aa_profile, dict_num_lineage = init_num_lineages('accession',
                                                                            'input/consensus_mutation_profile.tsv')
                    df_check_consensus_mutations = consensus_mutations(
                        lineage_list, df_dna_aa_profile)
                    df_check_consensus_mutations['used_mutation'] = df_created_consensus_mutations['dna_mutations']
                    df_check_consensus_mutations['check'] = df_check_consensus_mutations.apply(lambda row: compare_lists(
                        df_created_consensus_mutations[df_created_consensus_mutations['lineage'] == row['lineage']]['dna_mutations'].values[0], row['dna_mutations']), axis=1)
                    
                    print(df_check_consensus_mutations)
                    df_check_consensus_mutations.to_excel("output/consensus/mutation_consensus_check.xlsx")
                    '''
                    if True in df_check_consensus_mutations.check.values:
                        false_comparison = df_check_consensus_mutations.loc[df_check_consensus_mutations['column_name'] == True]
                        print("There are false mutations..")
                        print(false_comparison)
                    else:
                        print("No false mutations are found after Covsonar Check")
                    '''
                else:
                    print("No further check whether consensus have the right mutations")
                
            else:
                print("No option for Output type.")
                print("Mutation profile for nucleotides and amino acids: ", "\n", df_dna_aa_profile)
        else:
            print("Given path does not exist")

if __name__ == '__main__':
    main()
