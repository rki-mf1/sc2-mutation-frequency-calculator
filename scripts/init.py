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
import requests
import csv

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
        lineage_txt = mutation_file
        list_lineages_txt = Path(
            lineage_txt).read_text().replace('\n', '').split(",")
        pangolin_lineages = [i.split(",")[1] for i in Path(
            lineage_txt).read_text().split("\n")[:-1]]
        mutation_string = ' '.join(pangolin_lineages)
    if "aa" in mutation_file:
        # read .txt in list of strings
        # 1) aa_changes.txt
        # have to consider single/combinations of mutations
        aa_change_txt = mutation_file
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
    df_mutation_profile = pd.read_table(mutation_profile, low_memory=False)
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
            df_lin_frequency[lineage] = round(pd.DataFrame.from_dict(dict_num_lineage_mutations[lineage], orient='index') / num_lineages[lineage], 3)
            list_df_lin_frequency.append(df_lin_frequency)
    merged_df = reduce(lambda left, right: pd.merge(
        left, right, left_index=True, right_index=True, how='outer'), list_df_lin_frequency)
    return merged_df


########## optional parameters ############

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
                    consensus_seq_snp = consensus_seq
                    match = re.match(
                        r"([a-z]+)([0-9]+)([a-z]+)", mutation, re.I)
                    if match:
                        items = list(match.groups())
                        ref_position = int(items[1])
                        ref_base = consensus_seq[int(items[1])-1]
                        ref_mutation_base = items[0]
                        # warning message
                        if ref_base != ref_mutation_base:
                            print(
                                "Warning: reference nucleotide doesnt match the snp mutation")
                            print(f"Ref Seq at position {ref_position} is {ref_base}, ref base from mutation is {ref_mutation_base}")
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
                        if consenus_seq_insert[int(items[1])-1] != items[0]:
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
                {'lineage': lineage, 'dna_mutations': sorted(list_mutations)}, ignore_index=True)
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
    parser = argparse.ArgumentParser(description='SARS-CoV2 Mutation frequency calculator')
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
    parser.add_argument('-warning', '--mutation_warning', required=False, action='store_true',
                        help='optional: Warning message if first position of mutation is not the position in ref seq')
    # https://bioinformaticshome.com/db/collection/phylogenetics
    parser.add_argument('-vis', '--visualization', metavar='', required=False,
                        help='optional: visualization of consensus sequences')
    parser.add_argument('-cut', '--cut_off', metavar='', required=False, type=restricted_float,
                        help='optional: set threshold for aa mutations to create matrix or consensus')
    parser.add_argument('-upper', '--upper_cut', metavar='', required=False, type=restricted_float,
                        help='optional: set threshold for upper boarder to decide signature mutations')
    parser.add_argument('-lower', '--lower_cut', metavar='', required=False, type=restricted_float,
                        help='optional: set threshold for lower boarder to decide singature mutations')
    parser.add_argument('-sig', '--signature', metavar='', required=False,
                        help='optional: for given (set of) lineage return the set of signature mutations')
    parser.add_argument('-out', '--output', metavar='', required=False,
                        help='optional: set filename for Covsonar output')
    parser.add_argument('-level', '--mutation_level', metavar='', required=False,
                        help='optional: choose mutation level to get either aa or nt mutation frequencies')
    args = parser.parse_args()
    
    #load yml file 
    with open("config/config.yml", "r") as ymlfile:
        config = yaml.full_load(ymlfile)
    
    # Step1: check if covsonar mutation profile (output tsv) is given
    #AA Mutation als Fall für tsv generation DATAFRAME MERGE ÄNDERN AUF OUTTER
    if not args.mutation_profile:  
        print("Mutation profile from Covsonar doesnt exist:")
        if args.signature and ".xlsx" in args.signature:
            df_matrix = pd.read_excel(args.signature, index_col=0) 

            if args.lineages:
                file = open(args.lineages, "r")
                data = list(csv.reader(file, delimiter="\t"))
                file.close()
                lineages = [row[0].strip() for row in data]
            else:
                lineages = list(df_matrix.columns)
            
            if not args.upper_cut:  # cut-off row wise
                print("cut-off = 0.75 since no cut-off is selected")
                upper_cut = 0.75
            else:
                upper_cut = args.upper_cut
                print(f"cut-off is set to {upper_cut}")

            if not args.lower_cut:
                lower_cut = upper_cut/3
            else:
                lower_cut = args.lower_cut

            print("Signature mutations are: ")
            for lineage in lineages:
                candidates = df_matrix[lineage][(df_matrix[lineage] >= upper_cut)]
                candidate_mutations = list(candidates.index)
                df_matrix_candidates = df_matrix.loc[candidate_mutations].drop(lineage, axis=1)
                df_matrix_signature = df_matrix_candidates[(df_matrix_candidates <= lower_cut).all(1)]
                print(lineage, ':', list(df_matrix_signature.index))

        else:
            print("Covsonar running ...")
            
            if args.date_range:
                date_range = args.date_range
            else:
                print("A date range and database should be available for Covsonar to create mutation profile")

            if args.lineages and not args.aa_mutations:
                print("lineages but no aa")
                if isinstance(args.lineages, str):
                    lineages_string = args.lineages
                else:
                    lineages_string = txt_to_string(args.lineages)
                if args.output:
                    outfile = 'input/' + args.output
                else:
                    outfile = 'input/mutation_lineages_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match --lineage {lineages_string} --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}. \n"
                    "Choosen Lineages: \n"
                    f"{lineages_string}")
            elif args.aa_mutations and not args.lineages: #für jede mutation (single als auch comb) muss covsonar einzeln ausgeführt werden, freq berechnung für single/comb dann einzeln basierend auf merged tsv
                print("aa but no lin")
                aa_single_changes, aa_comb_changes = txt_to_string(args.aa_mutations)[0].replace(
                    " ", " " + '-i' + " "), txt_to_string(args.aa_mutations)[1]
                if args.output:
                    outfile = 'input/' + args.output
                else:
                    outfile = 'input/mutation_aa_single_profile.tsv'
                #comb_outfiles = 'input/combinations/mutation_aa{counter}_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match -i {aa_single_changes} --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}. \n"
                    "Choosen Amino acid Changes: \n"
                    f"{aa_single_changes}")
            elif args.lineages and args.aa_mutations:
                print("sowohl lin als auch aa")
                lineages_string = txt_to_string(args.lineages)
                aa_single_changes = txt_to_string(args.aa_mutations)[0].replace(
                    " ", " " + '-i' + " ")
                if args.output:
                    outfile = 'input/' + args.output
                else:
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
                    if args.output:
                        outfile = 'input/' + args.output
                    else:
                        outfile = 'input/mutation_profile.tsv'
                    os.system(
                        f"python3 covsonar/sonar.py match --date {args.date_range} --db {args.database} --tsv > {outfile}")
                    print(f"Covsonar created mutation profile in {outfile}.")
    
    else: # db and date not more relevant;
        if os.path.exists(args.mutation_profile):
            print("A mutation profile (covsonar output) exists")
            print("Read " + str(args.mutation_profile))

            df_mutation_profile = pd.read_table(args.mutation_profile, low_memory=False)
            date_range = str(df_mutation_profile.iloc[0]["date"] + ":" + df_mutation_profile.iloc[-1]["date"])

            df_dna_aa_profile, dict_num_lineage = init_num_lineages(
                'lineage', args.mutation_profile)
            
            print("Compute mutation frequency ...")
            
            # compute absolute values of mutations (dict)
            if args.lineages:  # .txt has to be string
                lineage_list = txt_to_string(args.lineages).split(' ')
            else:
                lineage_list = list(dict_num_lineage.keys())

            if args.matrix or args.signature: 
                print("Matrix will be created on the fly..")
                if args.mutation_level == "aa":
                    df_lineage_mutation_frequency = lineage_mutation_frequency(
                    "aa_profile", df_dna_aa_profile, lineage_list, dict_num_lineage)
                elif args.mutation_level == "nt":
                    df_lineage_mutation_frequency = lineage_mutation_frequency(
                        "dna_profile", df_dna_aa_profile, lineage_list, dict_num_lineage)
                
                if not args.cut_off: # rowwise
                    print("no cutoff")
                    if not args.aa_mutations:
                        print("and no aa mutations")
                        df_matrix = df_lineage_mutation_frequency
                    else:
                        print("but aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                else: #cutoff
                    print("cutoff is selected")
                    cut_off = args.cut_off
                    if not args.aa_mutations:
                        print("but no aa mutations")
                        df_matrix = df_lineage_mutation_frequency.fillna(
                            0)[(df_lineage_mutation_frequency.fillna(0) >= cut_off).any(axis=1)]
                    else: #cutoff + aa 
                        print("and aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix_aa = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                        df_matrix = df_matrix_aa.fillna(0)[(df_matrix_aa.fillna(0) >= args.cut_off).any(axis=1)]

                if args.matrix:
                    df_matrix = df_matrix.reindex(sorted(df_matrix.columns), axis=1)
                    df_matrix = df_matrix.iloc[df_matrix.index.map(lambda x: (x.split(':', 1)[0], int(
                        ''.join(filter(str.isdigit, x.split(':',1)[1]))))).argsort()]
                    matching_keys = set(dict_num_lineage.keys()).intersection(df_matrix.columns)
                    dict_filtered = {k: dict_num_lineage[k] for k in matching_keys}
                    dict_filtered_sorted = dict(sorted(dict_filtered.items()))
                    genome_number = pd.Series(
                        dict_filtered_sorted, name='Number of sequences detected')
                    df_matrix = pd.concat([genome_number.to_frame().T, df_matrix])
                    print(df_matrix)
                    
                    if args.output:
                        outfile = 'output/matrix/' + args.output
                    else:
                        outfile = 'output/matrix/frequency_matrix.xlsx'
                    df_matrix.to_excel(outfile)
                    print(f"Frequency matrix is created in {outfile}")

                elif args.signature:
                    if args.lineages:
                        file = open(args.signature, "r")
                        data = list(csv.reader(file, delimiter="\t"))
                        file.close()
                        sig_lineages = [row[0].strip() for row in data]
                    else:
                        lineages = list(df_matrix.columns)

                    if not args.upper_cut:  # cut-off row wise
                        print("cut-off = 0.75 since no cut-off is selected")
                        upper_cut = 0.75
                    else:
                        upper_cut = args.upper_cut
                        print(f"cut-off is set to {upper_cut}")

                    if not args.lower_cut:
                        lower_cut = upper_cut/3
                    else:
                        lower_cut = args.lower_cut

                    print("Signature mutations are: ")
                    for lineage in sig_lineages:
                        candidates = df_matrix[lineage][(df_matrix[lineage] >= upper_cut)]
                        candidate_mutations = list(candidates.index)
                        df_matrix_candidates = df_matrix.loc[candidate_mutations].drop(lineage, axis=1)
                        df_matrix_signature = df_matrix_candidates[(df_matrix_candidates <= lower_cut).all(1)]
                        print(lineage, ':', list(df_matrix_signature.index))

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
                    lineages_check = df_dna_aa_profile["accession"].tolist()
                    df_check_consensus_mutations = consensus_mutations(
                        lineage_list, df_dna_aa_profile)
                    df_check_consensus_mutations = df_check_consensus_mutations[df_check_consensus_mutations["lineage"].isin(lineages_check)].reset_index(drop=True)
                    df_check_consensus_mutations['used_mutation'] = df_created_consensus_mutations['dna_mutations']
                    df_check_consensus_mutations['check'] = df_check_consensus_mutations.apply(lambda row: compare_lists(
                        df_created_consensus_mutations[df_created_consensus_mutations['lineage'] == row['lineage']]['dna_mutations'].values[0], row['dna_mutations']), axis=1)
                    
                    print(df_check_consensus_mutations)

                    if False in df_check_consensus_mutations.check.values:
                        print("There are false mutations..")
                        false_comparison = df_check_consensus_mutations.loc[df_check_consensus_mutations['check'] == False].reset_index(
                            drop=True)
                        print(false_comparison)
                        false_comparison_lineage = false_comparison.loc[0]['lineage']
                        false_comparison_dna = sorted(false_comparison.loc[0]['dna_mutations'])
                        false_comparison_used = false_comparison.loc[0]['used_mutation']
                        print(
                            f"for lineage {false_comparison_lineage} among {false_comparison_dna} and {false_comparison_used}")
                        diff = set(false_comparison_dna) ^ set(false_comparison_used)
                        print(f"These differences are {diff}")
                    else:
                        print("Consensus sequences are Covsonar approved..")
                else:
                    print("No further check whether consensus have the right mutations")
            else:
                print("No option for Output type.")
                print("Mutation profile for nucleotides and amino acids: ", "\n", df_dna_aa_profile)
        else:
            print("Mutation profiles path does not exist")


if __name__ == '__main__':
    main()
