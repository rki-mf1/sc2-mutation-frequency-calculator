#!/usr/bin/env python3

# Master thesis FU Berlin
# RKI MF1
# Ashkan Ghassemi
# SARS-CoV-2 mutation frequency calculator - A covSonar utility tool

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
import xlsxwriter
from itertools import combinations
from collections import defaultdict
import matplotlib.pyplot as plt
#from matplotlib_venn import venn2_unweighted
import math
from pango_aliasor.aliasor import Aliasor

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
        mutation_string = Path(
            lineage_txt).read_text().split("\n")
    if "aa" in mutation_file:
        # read .txt in list of strings
        # 1) aa_changes.txt
        # have to consider single/combinations of mutations
        aa_change_txt = mutation_file
        list_aa_txt = Path(aa_change_txt).read_text().replace(
            '\n', ',').split(",")
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

def colToExcel(col):  # col is 1 based
    excelCol = str()
    div = col
    while div:
        (div, mod) = divmod(div-1, 26)  # will return (x, 0 .. 25)
        excelCol = chr(mod + 65) + excelCol

    return excelCol


def list_files_in_directory(path):
    # Check if the path exists
    if not os.path.exists(path):
        print(f"The path '{path}' does not exist.")
        return []

    # Use os.listdir to get a list of filenames in the directory
    filenames = os.listdir(path)

    # Filter out subdirectories and only keep files
    files = [filename for filename in filenames if os.path.isfile(
        os.path.join(path, filename))]

    return files
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
    list_df_num_mutations = []
    for lineage in lineages:
        df_lineage_subprofile = df_mutation_profile.loc[
            df_mutation_profile['lineage'] == lineage][mutation_type].reset_index(drop=True)
        df_lineage_subprofile = df_lineage_subprofile.dropna()
        result_list = [
            item for sublist in df_lineage_subprofile.str.split() for item in sublist]
        df_num_mutations = pd.Series(result_list).value_counts(
            ).rename_axis('mutation').reset_index(name='counts')
        df_num_mutations.set_index('mutation', inplace=True)
        df_num_mutations.index.name = None
        df_num_mutations.rename(columns={'counts': lineage}, inplace=True)
        df_num_mutations[lineage] = round(
            df_num_mutations[lineage] / num_lineages[lineage], 3)
        list_df_num_mutations.append(df_num_mutations)
    merged_df = reduce(lambda left, right: pd.merge(
        left, right, left_index=True, right_index=True, how='outer'), list_df_num_mutations)
    return merged_df

######### Matrix optimization ############

def orf_frameshift_matrix(matrix): #IFF there are duplicates in ORF1a and ORF1b
    '''Frameshift of ORF1a and ORF1b, both share the same mutation and have to be merged to ORF1ab
    input: df_matrix_noframeshift_unsorted
    output: df_matrix_frameshift_unsorted
    '''
    orf_matrix = matrix[matrix.index.str.startswith('ORF1a')
                        | matrix.index.str.startswith('ORF1b')].sort_index()
    index_names = orf_matrix.index.map(
        lambda x: x.split(':', 1)[1])
    duplicated_names = index_names[index_names.duplicated(keep=False)].unique()
    duplicated_df = orf_matrix[orf_matrix.index.str.contains(
        '|'.join(duplicated_names))]
    duplicate_indice_names = duplicated_df.index.tolist()
    duplicated_df_half = duplicated_df.iloc[:len(
        duplicated_df) // 2]
    duplicated_df_half.rename(index=lambda x: x.replace(
        x.split(':', 1)[0], 'ORF1ab', 1), inplace=True)
    df_matrix_noframeshift = matrix.drop(
        index=duplicate_indice_names)
    df_matrix = pd.concat(
        [df_matrix_noframeshift, duplicated_df_half])
    return df_matrix

def sort_matrix_columnandindex(matrix, num_lineages):
    '''Sort matrix according to columnnames (alphabetically) and indexnames by gene and position
    input: unsorted matrix, dict with num of lineages
    output: sorted matrix
    '''
    df_matrix = matrix.reindex(
        sorted(matrix.columns), axis=1)
    df_matrix = df_matrix.iloc[df_matrix.index.map(lambda x: (x.split(':', 1)[0], int(
        ''.join(filter(str.isdigit, x.split(':', 1)[1]))))).argsort()]
    #sample size
    matching_keys = set(num_lineages.keys()).intersection(
        df_matrix.columns)
    dict_filtered = {
        k: num_lineages[k] for k in matching_keys}
    dict_filtered_sorted = dict(sorted(dict_filtered.items()))
    genome_number = pd.Series(
        dict_filtered_sorted, name='Number of sequences detected')
    df_matrix = pd.concat(
        [genome_number.to_frame().T, df_matrix])
    return df_matrix

def sort_matrix_nt_columnandindex(matrix, num_lineages):
    '''Sort matrix according to columnnames (alphabetically) and indexnames by gene and position
    input: unsorted matrix, dict with num of lineages
    output: sorted matrix
    '''
    df_matrix = matrix
    df_matrix['Prefix'] = [
        re.match(r'([A-Za-z]+)', idx).group(1) for idx in df_matrix.index]
    df_matrix['Number'] = [int(re.search(r'(\d+)', idx).group(1))
                                               for idx in df_matrix.index]
    df_matrix = df_matrix.sort_values(
        by=['Prefix', 'Number'])
    df_matrix = df_matrix.drop(
        columns=['Prefix', 'Number'])
    
    #sample size
    matching_keys = set(num_lineages.keys()).intersection(
        df_matrix.columns)
    dict_filtered = {
        k: num_lineages[k] for k in matching_keys}
    dict_filtered_sorted = dict(sorted(dict_filtered.items()))
    genome_number = pd.Series(
        dict_filtered_sorted, name='Number of sequences detected')
    df_matrix = pd.concat(
        [genome_number.to_frame().T, df_matrix])
    return df_matrix

def init_num_labs(database, matrix, num_labs):
    '''Add the number of labs per lineage after the number of samples detected
    input: matrix and number of labs
    output: matrix with additional row at second position
    '''
    if database == "desh":
        matrix.loc['Labcounts'] = [num_labs.get(lineage, 0) for lineage in matrix.columns]
    else:
        matrix.loc['Countries'] = [num_labs.get(lineage, 0) for lineage in matrix.columns]
    labcounts_row = matrix.iloc[-1]
    matrix = matrix.drop(matrix.index[-1])
    sequences_index = matrix.index.get_loc(
        'Number of sequences detected')
    matrix = pd.concat([matrix.iloc[:sequences_index + 1],
        labcounts_row.to_frame().T, matrix.iloc[sequences_index + 1:]])
    return matrix

def pangolin_parent(input_string):
    recombinant = "X"
    if not input_string.startswith(recombinant):
        parts = input_string.split('.')

        if len(parts) > 4: #pangolin nomenclature 
            result = '.'.join(parts[:4])
        else:
            result = input_string
    else:  # https://github.com/MDU-PHL/pango-collapse/issues/7
        result = "XBB"
    return result

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

def create_consensus(infile, outdir, lineage_dna_frequency, threshold):
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

            file_out = f'{outdir}/sc2_consensus2_{lineage}.fasta'
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
                    consensus_seq, id=consensus_id, description="_Consensus_squared")

            with open(file_out, "w") as f:
                SeqIO.write(seq_record, f, "fasta")
            print(f"Consensus file for {consensus_id} is created in {file_out}")
    print("Consensus sequences can now be analysed using MSA and phylogenetic trees")
    return df_consensus_mutations

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

def bootstrap_value(reference, list_of_lists, iteration_size=100):
    ''' helperfunction for performing random bootstrap to calculate bootstrap-value per samplesize 
    input:
        reference list: ground truth mutation profile of every sample per lineage, 
        list_of_lists: list containing all mutation profiles from iteration step 
        iteration_size: how many mutation profiles to generate for each iteration sample size (default 100)  
    output: boostrap value 
    '''
    equal_count = 0
    for lst in list_of_lists:
        if compare_lists(reference, lst):
            equal_count += 1
    equal_ratio = equal_count / iteration_size
    return equal_ratio

def main():
    parser = argparse.ArgumentParser(description='SARS-CoV2 Mutation frequency calculator')
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
    parser.add_argument('-g', '--gradient', required=False, action='store_true',
                        help='optional: apply 2-color conditional formatting to create gradient for xlsx file')
    parser.add_argument('-lplot', '--lab_plot', required=False, action='store_true',
                        help='optional: plot about how many lineages coming from how many laboratories')
    parser.add_argument('-con', '--consensus', required=False, action='store_true',
                        help='optional: will create consensus sequences based on given lineages')
    parser.add_argument('-rb', '--bootstrap', required=False, action='store_true',
                        help='optional: will apply ranom bootstrap method to determine sufficient cut-off for minimum sample size per lineage')
    parser.add_argument('-tcd', '--treeshrink_cirlce_diagram', required=False, action='store_true',
                        help='optional: will plot a circle diagram, to compare phylotrees based on UShER and C^2 by checking which lineages are removed by treeshrink')
    parser.add_argument('-check', '--consensus_check', required=False, action='store_true',
                        help='optional: alterate path to check whether the created consensus seqs have the right mutations')
    parser.add_argument('-warning', '--mutation_warning', required=False, action='store_true',
                        help='optional: Warning message if first position of mutation is not the position in ref seq')
    parser.add_argument('-cut_freq', '--cut_off_frequency', metavar='', required=False, type=restricted_float,
                        help='optional: set threshold for aa mutations to create matrix or consensus')
    parser.add_argument('-cut_lin', '--cut_off_lineage', metavar='', required=False, type=int,
                        help='optional: set threshold for number of samples sequenced from a lineage')
    parser.add_argument('-sig', '--signature', required=False, action='store_true',
                        help='optional: for given (set of) lineage return the set of signature mutations')
    parser.add_argument('-is', '--ignore_sublineages', required=False, action='store_true',
                        help='optional: will ignore sublineages, to compute signature mutations for lineage-complexes')
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
        print("Mutation profile from Covsonar doesn't exist:")
        print("Covsonar running ...")
            
        if args.date_range:
            date_range = args.date_range
        else:
            print("Input missing: A date range and database should be available for Covsonar to create mutation profile!")

        if args.lineages and not args.aa_mutations:
            print("Input: lineages but no aa mutations")
            if isinstance(args.lineages, str):
                lineages_string = args.lineages
            else:
                lineages_string = txt_to_string(args.lineages)
            if args.output:
                outfile = 'input/tsv/' + args.output
            else:
                outfile = 'input/tsv/mutation_lineages_profile.tsv'
            os.system(
                f"python3 covsonar/sonar.py match --lineage {lineages_string} --date {args.date_range} --db {args.database} --tsv > {outfile}")
            print(f"Covsonar created mutation profile in {outfile}. \n"
                "Choosen Lineages: \n"
                f"{lineages_string}")
        elif args.aa_mutations and not args.lineages: #für jede mutation (single als auch comb) muss covsonar einzeln ausgeführt werden, freq berechnung für single/comb dann einzeln basierend auf merged tsv
            print("Input: aa mutations but no lineages")
            aa_single_changes, aa_comb_changes = txt_to_string(args.aa_mutations)[0].replace(
                " ", " " + '-i' + " "), txt_to_string(args.aa_mutations)[1]
            if args.output:
                outfile = 'input/tsv/' + args.output
            else:
                outfile = 'input/tsv/mutation_aa_single_profile.tsv'
            #comb_outfiles = 'input/combinations/mutation_aa{counter}_profile.tsv'
            os.system(
                f"python3 covsonar/sonar.py match -i {aa_single_changes} --date {args.date_range} --db {args.database} --tsv > {outfile}")
            print(f"Covsonar created mutation profile in {outfile}. \n"
                "Choosen Amino acid Changes: \n"
                f"{aa_single_changes}")
        elif args.lineages and args.aa_mutations:
            print("Input: lineages and aa mutations")
            lineages_string = txt_to_string(args.lineages)
            aa_single_changes = txt_to_string(args.aa_mutations)[0].replace(
                " ", " " + '-i' + " ")
            if args.output:
                outfile = 'input/tsv/' + args.output
            else:
                outfile = 'input/tsv/mutation_aa_lineages_profile.tsv'
            os.system(
                f"python3 covsonar/sonar.py match -i {aa_single_changes} --lineage {lineages_string} --date {args.date_range} --db {args.database} --tsv > {outfile}")
            print(f"Covsonar created mutation profile in {outfile}. \n"
                "Choosen Amino acid Changes: \n"
                f"{aa_single_changes} \n"
                "Choosen Lineages: \n"
                f"{lineages_string}")
        else:
            if not args.date_range or not args.database:
                print("Input missing: A date range and database should be available for Covsonar to create mutation profile")
            else:
                if args.output:
                    outfile = 'input/tsv/' + args.output
                else:
                    outfile = 'input/tsv/mutation_profile.tsv'
                os.system(
                    f"python3 covsonar/sonar.py match --date {args.date_range} --db {args.database} --tsv > {outfile}")
                print(f"Covsonar created mutation profile in {outfile}.")
    
    else: # db and date not more relevant;
        if os.path.exists(args.mutation_profile):
            print("A mutation profile (covsonar output) exists")
            print("Read " + str(args.mutation_profile))

            df_mutation_profile = pd.read_table(args.mutation_profile, low_memory=False)
            df_mutation_profile['date'] = pd.to_datetime(df_mutation_profile['date'])
            min_date, max_date = df_mutation_profile['date'].min(), df_mutation_profile['date'].max()
            date_range = f"{min_date.strftime('%Y-%m-%d')}:{max_date.strftime('%Y-%m-%d')}"
            print("Time Period:", date_range)
                        
            print("Absolut number of sequenced lineages will be counted...")
            df_dna_aa_profile, dict_num_lineage = init_num_lineages('lineage', args.mutation_profile)

            if args.cut_off_lineage:
                sample_cut_off = args.cut_off_lineage
                dict_filter_num_lineage = {k: v for k, v in dict_num_lineage.items() if v >= args.cut_off_lineage}
                dict_other_lineages = {
                    k: v for k, v in dict_num_lineage.items() if v < args.cut_off_lineage} 
            else:
                sample_cut_off = 10
                dict_filter_num_lineage = {k: v for k, v in dict_num_lineage.items() if v >= sample_cut_off}
                dict_other_lineages = {
                    k: v for k, v in dict_num_lineage.items() if v < sample_cut_off}
            
            if args.lineages:  # .txt has to be string
                tmp = list(dict_filter_num_lineage.keys())
                input_list = txt_to_string(args.lineages)
                lineage_list = list(set(input_list) & set(tmp))
            else:
                lineage_list = list(dict_filter_num_lineage.keys())
            
            sorted_lineage_list = sorted(lineage_list)
            dict_filter_num_lineage = {key: dict_filter_num_lineage[key] for key in lineage_list}
            dict_filter_num_lineage = dict(sorted(dict_filter_num_lineage.items()))
            #print("Number of samples per lineage:", dict_filter_num_lineage)
            #print(dict_other_lineages)
            
            print("Compute mutation frequency ...") 

            if args.matrix or args.signature: 
                print("Matrix will be created on the fly..")

                # column "other lineages"
                lineage_list.append("other_lineages")
                other_lineages_samples = sum(dict_other_lineages.values())
                dict_filter_num_lineage['other_lineages'] = other_lineages_samples
                df_dna_aa_profile['lineage'] = df_dna_aa_profile['lineage'].apply(
                    lambda x: 'other_lineages' if x in dict_other_lineages else x)
                df_mutation_profile['lineage'] = df_mutation_profile['lineage'].apply(
                    lambda x: 'other_lineages' if x in dict_other_lineages else x)
                
                #labdiversity
                #country info is in column "zip"
                
                database = "gisaid"
                if not df_mutation_profile['lab'].isnull().values.all():
                    database = "desh"
                    selected_mutation_profile = df_mutation_profile.loc[df_mutation_profile['lineage'].isin(
                        lineage_list), ['lab', 'lineage']]
                    lab_zip_counts = selected_mutation_profile.groupby('lineage')[
                    'lab'].nunique()
                    
                    for lineage, lab_count in lab_zip_counts.items():
                        if lab_count <= 1:
                            print(
                                f"Warning: {lineage} comes from only {lab_count} labs! Note, that {lineage} is not excluded from further analysis but should be considered critically.")
                else:
                    selected_mutation_profile = df_mutation_profile.loc[df_mutation_profile['lineage'].isin(
                        lineage_list), ['zip', 'lineage']]
                    lab_zip_counts = selected_mutation_profile.groupby('lineage')[
                    'zip'].nunique()

                if args.mutation_level == "aa":
                    df_lineage_mutation_frequency = lineage_mutation_frequency(
                    "aa_profile", df_dna_aa_profile, lineage_list, dict_filter_num_lineage)
                elif args.mutation_level == "nt":
                    df_lineage_mutation_frequency = lineage_mutation_frequency(
                        "dna_profile", df_dna_aa_profile, lineage_list, dict_filter_num_lineage)
                else:
                    "MissingError: mutation level is not defined"
                
                if not args.cut_off_frequency: # rowwise
                    print("no cutoff")
                    if not args.aa_mutations:
                        print("and no aa mutations")
                        df_matrix_noframeshift_unsorted = df_lineage_mutation_frequency
                    else:
                        print("but aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix_noframeshift_unsorted = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                
                else: #cutoff
                    print("cutoff is selected")
                    cut_off_frequency = args.cut_off_frequency
                    if not args.aa_mutations:
                        print("but no aa mutations")
                        df_matrix_noframeshift_unsorted = df_lineage_mutation_frequency.fillna(
                            0)[(df_lineage_mutation_frequency.fillna(0) >= cut_off_frequency).any(axis=1)]
                    else: #cutoff + aa 
                        print("and aa mutations")
                        list_aa_single_changes = txt_to_string(
                            args.aa_mutations)[0].split(" ")
                        df_matrix_aa = df_lineage_mutation_frequency[df_lineage_mutation_frequency.index.isin(
                            list_aa_single_changes)]
                        df_matrix_noframeshift_unsorted = df_matrix_aa.fillna(0)[(df_matrix_aa.fillna(0) >= args.cut_off_frequency).any(axis=1)]
                
                df_matrix_frameshift_unsorted = orf_frameshift_matrix(
                    df_matrix_noframeshift_unsorted)
                
                if args.mutation_level == "aa":
                    df_matrix_frameshift_unsorted_gene_order = sort_matrix_columnandindex(
                        df_matrix_frameshift_unsorted, dict_filter_num_lineage)
                    # also sort mutations according to gene presence in the genome (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8067447/, Fig.1)
                    
                    gene_order = ['ORF1a', 'ORF1ab','ORF1b', 'S', 'ORF3a', 'ORF3d', 'E', 'M',
                                  'ORF6', 'ORF7a', 'ORF7b', 'ORF8', 'ORF9b', 'ORF14', 'N', 'ORF10']
                    aa_mutations = df_matrix_frameshift_unsorted_gene_order.index.to_list()[1:]
                    
                    def sort_gene_mutations(item):
                        gene_name = item.split(':')[0]  # Extract the gene name from the item
                        try:
                            # Return the index of the gene name in the desired order
                            return gene_order.index(gene_name)
                        except ValueError:
                            return len(gene_order)
                    
                    sorted_mutations = ['Number of sequences detected'] + sorted(aa_mutations, key=sort_gene_mutations)
                    df_matrix_without_lab = df_matrix_frameshift_unsorted_gene_order.reindex(
                        sorted_mutations)
                    
                elif args.mutation_level == "nt":
                    df_matrix_without_lab = sort_matrix_nt_columnandindex(
                        df_matrix_frameshift_unsorted, dict_filter_num_lineage)
                
                df_matrix_without_parent = init_num_labs(database, df_matrix_without_lab, lab_zip_counts)
                
                #parent-child relationship
                sublineages_list = df_matrix_without_parent.columns.to_list()
                
                aliasor = Aliasor() #https: // github.com/corneliusroemer/pango_aliasor
                aliasor_uncompressed = [aliasor.uncompress(
                    sublineage) for sublineage in sublineages_list]
                parent_lineages = [pangolin_parent(lineage_uncompressed)
                                for lineage_uncompressed in aliasor_uncompressed]
                
                parent_lineages[parent_lineages.index('other_lineages')] = '---'
                parent_lineages_row = pd.DataFrame(
                    [parent_lineages], columns=df_matrix_without_parent.columns)
                parent_lineages_row.index = ['Parent lineage']
                df_matrix = pd.concat(
                    [parent_lineages_row, df_matrix_without_parent])
                
                if args.matrix:
                    print(df_matrix)
                    if args.output:
                        outfile = args.output
                    else:
                        outfile = 'output/matrix/' + \
                            date_range.replace(":", "_") + '.xlsx'
                    df_matrix.to_excel(outfile)
                    print(f"Frequency matrix is created in {outfile}")

                    if args.gradient:
                        excel_file = outfile
                        sheet_name = 'Mutation Frequency Table'
                        writer = pd.ExcelWriter(excel_file, engine='xlsxwriter')
                        df_matrix.to_excel(writer, sheet_name=sheet_name)
                        workbook = writer.book
                        worksheet = writer.sheets[sheet_name]
                        end = colToExcel(len(df_matrix.columns)+1) 
                        coordinate = 'B5:' + end + str(len(df_matrix)+1)
                        worksheet.conditional_format(coordinate, {'type': '2_color_scale',
                                             'min_color': 'white',
                                             'max_color': 'red'})
                        worksheet.write('A1', date_range)
                        writer.save()

                    if args.lab_plot:
                        # plot lab counts -> how many lineages come from 1, 2, 3, ... n labs
                        count_dict = defaultdict(int)
                        for value in lab_zip_counts:
                            count_dict[value] += 1
                        count_df = pd.DataFrame(
                            count_dict.items(), columns=['Value', 'Count'])
                        sorted_count_df = count_df.sort_values(
                            by='Value').reset_index(drop=True)

                        x_values, y_values = sorted_count_df['Value'], sorted_count_df['Count']
                        plt.bar(x_values, y_values)
                        plt.ylabel('Number of Lineages')
                        if database == "desh":
                            plt.xlabel('Number of Labs')
                            plt.title('Lab diversity')
                        else: 
                            plt.xlabel('Number of Countries')
                            plt.title('Country diversity')
                        
                        # deletion average frequencies against SNPs for all lineages
                    '''
                    if args.del_freq_plot:
                        df_matrix['Freq_Mean'] = df_matrix.mean()
                        print(df_matrix)
                        quit()
                        x_values, y_values = sorted_count_df['Value'], sorted_count_df['Count']
                        plt.bar(x_values, y_values)
                        plt.xlabel('Number of Labs')
                        plt.ylabel('Number of Lineages')
                        plt.title('Lab diversity')

                        dict_mut_freq = {}
                        for index, row in df_matrix.iloc[3:].iterrows():
                            row_values = row.tolist()
                            dict_mut_freq[index] = mean(row_values)

                        dict_del_mutations = {
                            key: value for key, value in dict_mut_freq.items() if "del" in key}
                        dict_snp_mutations = {
                            key: value for key, value in dict_mut_freq.items() if not "del" in key}

                        plt.bar(list(dict_del_mutations.keys()),
                                list(dict_del_mutations.values()))
                        plt.xlabel('Deletion mutations')
                        plt.ylabel('Frequency')
                        plt.title(
                            f'Frequency of deletions within {date_range}')
                        plt.xticks(rotation=90)
                        plt.show()

                        plt.bar(list(dict_snp_mutations.keys()),
                                list(dict_snp_mutations.values()))
                        plt.xlabel('Deletion mutations')
                        plt.ylabel('Frequency')
                        plt.title(
                            f'Frequency of deletions within {date_range}')
                        plt.xticks(rotation=90)
                        plt.show()

                        if args.output:
                            plt.savefig(f"output/matrix/{args.output}.png")
                        else:
                            plt.savefig(f"output/matrix/{date_range}_labdiversity.png")
                    '''
                elif args.signature:
                    # try out for input/tsv/2023-10-11_gisaid_ww-check_2023-09-11_2023-10-11.tsv
                    df_matrix = df_matrix.tail(-3)
                    lineages_of_interest = df_matrix.columns.to_list()

                    #if args.ignore_sublineages
                    #changes shape of matrix, less columns since sublineages will be removed

                    print(df_matrix)
                    columns_to_delete = ['HK.3', 'JG.3']
                    df_matrix.drop(columns=columns_to_delete, inplace=True)

                    #count values (in how many lineages a mutations is characteristic)
                    counts = (df_matrix >= args.cut_off_frequency).sum(axis=1)
                    frequencies = df_matrix[df_matrix >= args.cut_off_frequency].apply(
                        lambda row: row[row.notnull()].tolist(), axis=1)
                    #count and frequencies depict df_matrix
                    count_freq_table = pd.DataFrame(
                        {'count': counts, 'frequency': frequencies})
                    count_freq_table['frequency'] = count_freq_table['frequency'].apply(tuple)
                    count_freq_table['lineages'] = df_matrix.apply(lambda x: x.index[x >= args.cut_off_frequency].to_list(), axis=1)
                    count_freq_table_sort = count_freq_table.sort_values(by=['count', 'frequency'], ascending=[True, False])
                    count_freq_table_sort['frequency'] = count_freq_table_sort['frequency'].apply(list)

                    #count = 1
                    single_signature_table = count_freq_table_sort[count_freq_table_sort['count'] == 1]
                    single_signature_table['lineages'] = single_signature_table['lineages'].apply(lambda x: x[0])
                    unique_lineages = single_signature_table['lineages'].unique() #vorkommen können aus count tabelle raus
                
                    print("\n", "Unique lineages, die single signature haben: ", unique_lineages, "\n")
                    
                    print(single_signature_table)
                    df_signature_mutations = pd.DataFrame(
                        {'Lineages': single_signature_table['lineages'].unique()}).sort_values(by=['Lineages'])
                    df_signature_mutations['Combination Length'] = 1
                    df_signature_mutations['Signature Mutations'] = single_signature_table.groupby('lineages').apply(lambda x: x.index.tolist())
                    df_signature_mutations.reset_index(drop=True, inplace=True)
                    print(df_signature_mutations)
                    print(single_signature_table.groupby(
                        'lineages').apply(lambda x: x.index.tolist()))
                    
                    single_signatures = {} #aufheben, um mit combinations zu mergen und dann als .tsv/.txt speichern
                    for lineage in unique_lineages:
                        grouped = single_signature_table[single_signature_table['lineages'] == lineage].index.tolist()
                        if 'lineages' in single_signatures:
                            single_signatures['lineages'][lineage] = grouped
                        else:
                            single_signatures['lineages'] = {lineage: grouped}
                    
                    for column, value_grouped in single_signatures.items():
                        print(f"Single signature mutations for lineages within {date_range}")
                        for value, grouped in value_grouped.items():
                            print(f"Lineage: {value}")
                            print(grouped)
                    
                    #count > 1
                    
                    print(
                        "\n", "Unique lineages, die single signature haben: ", unique_lineages, "\n")
                    loi_find_combinations = [item for item in lineages_of_interest if item not in unique_lineages]
                    print("Lineages which need a combination of signature mutations: ",loi_find_combinations)
                    
                    combination_signature_table = count_freq_table_sort[count_freq_table_sort['count'] != 1]
                    combination_signature_filtered = combination_signature_table[~combination_signature_table['lineages'].apply(lambda x: set(x).issubset(unique_lineages))]
                    combination_signature_filtered[['frequency', 'lineages']] = combination_signature_filtered.apply(lambda row: pd.Series(
                        [list(frequency) for frequency in zip(*sorted(zip(row['frequency'], row['lineages']), reverse=True))]),
                        axis=1)
                    combination_signature_filtered['frequency'] = combination_signature_filtered['frequency'].apply(tuple)
                    combination_signature_filtered = combination_signature_filtered.sort_values(
                        by=['count', 'frequency'], ascending=[True, False])
                    combination_signature_filtered['frequency'] = combination_signature_filtered['frequency'].apply(list)
                    print(combination_signature_filtered)
                    
                    for loi in loi_find_combinations:
                        print(loi)
                        signature_combinations_loi = {}
                        df_loi_selected = combination_signature_filtered[combination_signature_filtered['lineages'].apply(
                            lambda x: loi in x)]
                        if not df_loi_selected.empty: #else: f"No mutation occured significantly in lineage: {loi}"
                            k = 2
                            condition_met = False
                            while k < 6:
                                k_combinations_loi = list(combinations(df_loi_selected.index, k))
                                for partition in k_combinations_loi:
                                    lineages = [df_loi_selected.loc[x, 'lineages'] for x in partition]
                                    intersection_all = list(set.intersection(*map(set, lineages)))
                                    if len(intersection_all) == 1 and intersection_all[0] == loi:
                                        print(f"{k}-combination: ")
                                        print(partition)
                                        condition_met = True
                                if condition_met:
                                    break
                                k+=1
                        
            
                    """
                        k = 2

                        if not df_loi_selected.empty:
                            pairwise_combinations_loi = list(combinations(
                                df_loi_selected.index, 2)) 
                    quit()
                    for loi in loi_find_combinations:    
                        df_loi_selected = combination_signature_filtered[combination_signature_filtered['lineages'].apply(
                            lambda x: loi in x)]
                        
                        if not df_loi_selected.empty:
                            pairwise_combinations_loi = list(combinations(
                                df_loi_selected.index, 2)) #n(n-1)/2
                            
                            signature_combinations_loi = []
                            for pair in pairwise_combinations_loi:
                                list1 = df_loi_selected.loc[pair[0], 'lineages']
                                list2 = df_loi_selected.loc[pair[1], 'lineages']
                                intersection = set(list1) & set(list2)
                                if len(intersection) == 1 and loi in intersection:
                                    signature_combinations_loi.append(pair)
                            if signature_combinations_loi:
                                print(loi)
                                print(signature_combinations_loi)
                            else:
                                print(
                                    f"no signature mutations for {loi} found")
                        else:
                            print(f"no signature mutations for {loi} found")
                    """

            elif args.consensus: 
                print("Mutation Frequency per lineage will be computed..")
                
                df_lineage_mutation_frequency = lineage_mutation_frequency(
                    "dna_profile", df_dna_aa_profile, sorted_lineage_list, dict_filter_num_lineage)
                num_of_lineages = len(sorted_lineage_list)
                
                if args.bootstrap:
                    print("random bootstrap method determines minimum sample size per lineage...")

                    bootstrap_threshold = 0.75
                    sample_size_factor = 1/2  
                    print(f"Bootstrap treshold: {bootstrap_threshold}, the sample size will be reduced by: {sample_size_factor}")

                    #dict_filter_num_lineage = dict(sorted(dict_filter_num_lineage.items(), key=lambda x: x[1], reverse=True))
                    #sorted_lineage_list = list(dict_filter_num_lineage.keys())
                    
                    df_bootstrap_lineages = pd.DataFrame(columns=['Lineages', 'Sample_size', 'Reduction_factor', 'Bootstrap_value'])
                    num_reflect_ground_truth = []
                    rat_reflect_ground_truth = []
                    for lineage in sorted_lineage_list:
                        print(f"---{lineage}---")
                        sample_size = dict_num_lineage[lineage]
                        df_lineage_mutation_profile = df_dna_aa_profile[
                            df_dna_aa_profile['lineage'] == lineage][['lineage', 'dna_profile']]
                        ground_truth = df_lineage_mutation_frequency.index[df_lineage_mutation_frequency[lineage]
                                                            >= bootstrap_threshold].tolist()
                        sample_size = math.floor(sample_size*sample_size_factor)
                        ground_truth = sorted(
                            ground_truth, key=lambda x: int(re.search(r'\d+', x).group()))
                        #print(ground_truth)
                    # presence/absence vector matrix for each sample of a lineage with aggregated count
                        #print(f"Ground truth for lineage {lineage}: ", ground_truth)
                        mutation_profile_aggregated = df_lineage_mutation_profile['dna_profile'].str.split(
                        )
                        df_presence_absence = pd.DataFrame([[col in row for col in ground_truth]
                                                            for row in mutation_profile_aggregated], columns=ground_truth, dtype=int)
                        df_presence_absence['Count'] = df_presence_absence.sum(
                            axis=1)
                        df_presence_absence = df_presence_absence.sort_values(
                            by='Count', ascending=False).reset_index(drop=True)
                        count_equal_length = (
                            df_presence_absence['Count'] == len(ground_truth)).sum()
                        print(df_presence_absence)
                        date_range = date_range.replace(":", "_")
                        if not os.path.exists(f"output/presence_absence/{date_range}/"):
                            os.makedirs(
                                f"output/presence_absence/{date_range}/")
                        df_presence_absence.to_csv(f"output/presence_absence/{date_range}/presence-absence_{lineage}_{date_range}.csv",index=False)
                        num_reflect_ground_truth.append(count_equal_length)
                        rat_reflect_ground_truth.append(
                            count_equal_length/dict_filter_num_lineage[lineage])
                        #print(f"Number of samples which exactly reflect the ground truth of {lineage}:", count_equal_length)
                        #print(f"Ratio of samples which exactly reflect the ground truth of {lineage}:", count_equal_length/dict_filter_num_lineage[lineage])

                    ''' #bootstrap approach
                    x_values, y_values = sorted_lineage_list, num_reflect_ground_truth
                    plt.bar(x_values, y_values)
                    plt.xlabel('Lineages')
                    plt.ylabel(
                        'Number of samples exactly reflect the ground truth')
                    plt.title('Determining minimum sample size')
                    overall_sample_size = list(
                        dict_filter_num_lineage.values())
                    for x, y, add_val in zip(x_values, y_values, overall_sample_size):
                        plt.text(x, y, str(add_val), ha='center', va='bottom')
                    plt.xticks(rotation=90)
                    plt.show()
                    a_values, b_values = sorted_lineage_list, rat_reflect_ground_truth
                    plt.bar(a_values, b_values)
                    plt.xlabel('Lineages')
                    plt.ylabel(
                        'Ratio of samples exactly reflect the ground truth over overall sample size')
                    plt.title('Determining minimum sample size')
                    overall_sample_size = list(
                        dict_filter_num_lineage.values())
                    for x, y, add_val in zip(a_values, b_values, overall_sample_size):
                        plt.text(x, y, str(add_val), ha='center', va='bottom')
                    plt.xticks(rotation=90)
                    plt.show()

                    
                        # random bootstrap approach
                        sample_size = math.floor(
                            sample_size*sample_size_factor)
                        dict_filter_num_lineage[lineage] = sample_size
                        
                        df_lineage_sample = pd.DataFrame(columns=['Lineages', 'Sample_size', 'Reduction_factor', 'Bootstrap_value'])
                        while sample_size >= 2:
                            bootstrap_range = 100
                            mutation_profile_sampled = []
                            for k in range(bootstrap_range):
                                df_lineage_mutation_profile_subsample = df_lineage_mutation_profile.sample(n=sample_size)['dna_profile'].reset_index(drop=True)
                                df_lineage_subprofile = df_lineage_mutation_profile_subsample.dropna()
                                mutation_profile_aggregated = [
                                    item for sublist in df_lineage_subprofile.str.split() for item in sublist]
                                df_num_mutations = pd.Series(mutation_profile_aggregated).value_counts(
                                ).rename_axis('mutation').reset_index(name='counts')
                                df_num_mutations.set_index('mutation', inplace=True)
                                df_num_mutations.index.name = None
                                df_num_mutations.rename(columns={'counts': lineage}, inplace=True)
                                df_num_mutations[lineage] = round(df_num_mutations[lineage] / dict_filter_num_lineage[lineage], 3)
                                list_mutations = df_num_mutations.index[df_num_mutations[lineage]
                                                                                     >= bootstrap_threshold].tolist()
                                mutation_profile_sampled.append(list_mutations)
                            
                            boostrap = bootstrap_value(ground_truth, mutation_profile_sampled, bootstrap_range)
                            #print(f"Bootstrap-value for {sample_size}: ", boostrap)
                            df_lineage_sample.loc[len(df_lineage_sample)] = [
                                lineage, sample_size, sample_size_factor, boostrap]
                            sample_size = math.floor(sample_size*sample_size_factor)
                            dict_filter_num_lineage[lineage] = sample_size
                        if (df_lineage_sample['Bootstrap_value'] >= 0.75).any():
                            row_minimum_sample_size = df_lineage_sample[df_lineage_sample['Bootstrap_value']>= 0.75].iloc[-1]
                            df_bootstrap_lineages = df_bootstrap_lineages.append(
                                row_minimum_sample_size, ignore_index=True)
                        else:
                            df_bootstrap_lineages.loc[len(df_bootstrap_lineages)] = [
                                lineage, df_lineage_sample['Sample_size'].iloc[0]*2, sample_size_factor, 1.0]
                    print(df_bootstrap_lineages)

                    date_range = date_range.replace(":", "_")
                    df_bootstrap_lineages.to_csv(f'output/bootstrap/{date_range}/2023-09-06_accession-all_min-sample-size_cut-off-0-75_{date_range}.tsv', sep="\t", index=False)
                    
                    #bootstrap plot
                    x_values, y_values = sorted_lineage_list, df_bootstrap_lineages['Sample_size']
                    plt.bar(x_values, y_values)
                    plt.xlabel('Lineages')
                    plt.ylabel('Minimum sample size')
                    plt.title('Random Bootstrap for determining minimum sample size') 
                    plt.savefig(
                        f"output/bootstrap/{date_range}/2023-09-06_accession-all_min-sample-size_cut-off-0-75.png")
                '''
                
                elif args.treeshrink_cirlce_diagram: 
                    directory_path = "output/consensus/2023-09-02_desh-accesions-all_cut-lin-5_consensus-squared/"
                    file_list = list_files_in_directory(directory_path)
                    if file_list:
                        list_all_intersected_lineages = [filename.split("_")[2].rsplit(".", 1)[
                            0] for filename in file_list]
                    else:
                        print("No files found in the directory.")
                    
                    df_treeshrink_removed_c1 = pd.read_table('input/phylotree_usher_c2/treeshrink_removed_c1.tsv', low_memory=False)
                    list_treeshrink_removed_c1 = sorted(df_treeshrink_removed_c1['lineage'].unique())
                    list_treeshrink_removed_c1_same = [
                        item for item in list_treeshrink_removed_c1 if item in list_all_intersected_lineages]
                    
                    file_path = "input/phylotree_usher_c2/removed_cut-lin-10.txt"
                    with open(file_path, 'r') as file:
                        # 2. Read the file's contents
                        file_contents = file.read()
                    list_treeshrink_removed_c2 = file_contents.split('\t')
                    if list_treeshrink_removed_c2[-1] == '\n':
                        list_treeshrink_removed_c2.pop()
                    list_treeshrink_removed_c2 = sorted(list_treeshrink_removed_c2)
                    list_treeshrink_removed_c2_same = [
                        item for item in list_treeshrink_removed_c2 if item in list_all_intersected_lineages]
                    
                    '''
                    venn2_unweighted([set(list_treeshrink_removed_c1), set(list_treeshrink_removed_c2)], ('C1', 'C2'))
                    plt.title(
                        "Removed lineages by treeshrink from UShER and C^2 phylotree")
                    plt.savefig("input/phylotree_usher_c2/venn_treeshrink-result.png")
                    
                    venn2_unweighted([set(list_treeshrink_removed_c1_same), set(list_treeshrink_removed_c2_same)], ('C1', 'C2'))
                    plt.title("Removed lineages by treeshrink from UShER and C^2 phylotree (same lineages)")
                    plt.savefig("input/phylotree_usher_c2/venn_treeshrink-result_same-lineages.png")
                    '''

                else: 
                    print("creating consensus sequences...")
                    file_in = 'data/NC_045512.2.fasta.txt'
                    
                    if not args.cut_off_frequency:  # cut-off columnwise
                        print("cut-off = 0.75 since no cut-off is selected")
                        cut_off_frequency = 0.75
                    else:
                        cut_off_frequency = args.cut_off_frequency
                        print(f"cut-off is set to {cut_off_frequency}")

                    fasta_dir = "output/consensus/"
                    if args.output:
                        outpath = fasta_dir + args.output
                    else:
                        date_range = date_range.replace(":", "_")
                        outpath = fasta_dir + date_range

                    if not os.path.exists(outpath + '/'):
                        os.makedirs(outpath + '/')

                    df_created_consensus_mutations = create_consensus(file_in, outpath, df_lineage_mutation_frequency, cut_off_frequency)
                    
                    records = []
                    for filename in os.listdir(outpath):
                        if filename.endswith(".fasta") or filename.endswith(".fa"):
                            filepath = os.path.join(fasta_dir, filename)
                            for record in SeqIO.parse(filepath, "fasta"):
                                records.append(record)

                    # Write the merged fasta sequences to a new file
                    if args.output:
                        multi_fasta = f"{args.output}/{args.output}_multi.fasta"
                    else:
                        date_range = date_range.replace(":", "_")
                        multi_fasta = f"{outpath}/{date_range}_multi.fasta"

                    with open(multi_fasta, "w") as f:
                        SeqIO.write(records, f, "fasta")
                    print(f"Merged consensus file is created in {multi_fasta} and contains {num_of_lineages}")
                    
                    if args.consensus_check:
                        print("consensus seq check if the mutations are right")
                        print("Created Consenus will be checked via Covsonar..")
                        print("Mutation profiles will be generated")

                        os.system('touch data/consensus.db')
                        os.system(
                            f'python3 covsonar/sonar.py add -f {outpath} --db data/consensus.db --cpus 8 --noprogress')
                        os.system(
                            'python3 covsonar/sonar.py match --db data/consensus.db --tsv > input/consensus_{outpath}.tsv')
                        
                        print("Mutations will be checked...")
                        df_dna_aa_profile, dict_filter_num_lineage = init_num_lineages('accession',
                                                                                'input/consensus_{outpath}.tsv')
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
