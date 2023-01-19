#!/usr/bin/env python
# coding: utf-8



import re
import os
import sys
import argparse
import pandas as pd 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO
from functions import  remove_short_sequeces,remove_redundant_sequeces,blast,generate_blast_hit_list,separate_hit_sequences,removal_of_sequences,blast_host



parser = argparse.ArgumentParser(prog="drgTFinder.py", formatter_class=argparse.RawTextHelpFormatter, description="""
drgTFinder is an automated pipeline for detecting Novel drug targets using proteome sequence from \n""",epilog="""
    Examples:
    python3 drgTFinder.py -db database -s test.fasta -o output -t 4\n""")
parser.add_argument('--sequence', '-s', type=str, default="", required=True, help='Whole proteome sequence file in FASTA format.')
parser.add_argument('--length', '-l', type=str, default="100", required=False, help='Sequence Lenth cutoff (default = 100)')
parser.add_argument('--cutoff', '-c', type=str, default="0.8", required=False, help='CD-Hit cutoff value (default = 0.8)')
parser.add_argument('--output', '-o', type=str, default="", required=False, help='Output directory')
parser.add_argument('--database', '-db', type=str, default="", required=False, help='Database directory')
parser.add_argument('--threads', '-t', type=str, default="4", required=False, help='Number of threads')









args = parser.parse_args()
file_name =  args.sequence
file_title = file_name.split(".")
if file_title[1] not in ["fasta","fa", "faa"] :
	sys.exit('Please provide a FASTA file')

file_input = file_name


# if args.output == "":
# 	sys.exit('please provide a output directory')


def main():

    os.system('mkdir temp > /dev/null')
    os.system("mkdir {}".format(args.output))

    remove_short_sequeces(args.length, file_input, "temp/large_seq.fasta")

    remove_redundant_sequeces("temp/large_seq.fasta","temp/non_redundant.fasta", args.cutoff)



    #non-host sequence
    blast_host("temp/non_redundant.fasta", "temp/blast_host.csv", args.database + "/host", args.threads )
    generate_blast_hit_list("temp/blast_host.csv", "temp/host_homologous.list")
    separate_hit_sequences(args.sequence, "temp/host_homologous.list", "temp/host_homologous.fasta")
    removal_of_sequences("temp/non_redundant.fasta", "temp/host_homologous.fasta", "temp/host_non_homologous.fasta")


    #essential sequence
    blast("temp/host_non_homologous.fasta", "temp/blast_deg.csv", args.database + "/deg", args.threads )
    generate_blast_hit_list("temp/blast_deg.csv", "temp/deg_homologous.list")
    separate_hit_sequences(args.sequence, "temp/deg_homologous.list", "temp/deg_homologous.fasta")


    #non-related-metabolic sequences
    blast("temp/deg_homologous.fasta", "temp/blast_kegg.csv", args.database + "/kegg", args.threads )
    generate_blast_hit_list("temp/blast_kegg.csv", "temp/kegg_homologous.list")
    separate_hit_sequences(args.sequence, "temp/kegg_homologous.list", "temp/kegg_homologous.fasta")



    os.system("cp temp/kegg_homologous.fasta {}/targets.fasta".format(args.output))

    #final_report(file_input)


#Go pipeline, go!!!
if __name__ == "__main__":
   main()


