#!/usr/bin/env python
# coding: utf-8

# In[1]:

import re
import os
import sys
import argparse
import pandas as pd 
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO



parser = argparse.ArgumentParser(prog="drgTFinder.py", formatter_class=argparse.RawTextHelpFormatter, description="""
drgTFinder is an automated pipeline for detecting Novel drug targets using proteome sequence from \n""",epilog="""
    Examples:
    python3 drgTFinder.py -s whole_proteome.fasta\n""")
parser.add_argument('--sequence', '-s', type=str, default="", required=True, help='Whole proteome sequence file in FASTA format.')
parser.add_argument('--type', '-t', type=str, default="", required=True, help='Bacteria type(p or n')
parser.add_argument('--length', '-l', type=str, default="100", required=False, help='Sequence Lenth cutoff (default = 100)')
parser.add_argument('--cutoff', '-c', type=str, default="0.8", required=False, help='CD-Hit cutoff value (default = 0.8)')
parser.add_argument('--localization', '-loc', type=str, default="CytoplasmicMembrane", required=False, help='Subcellular localization (Cytoplasmic, Extracellular, Cellwall, CytoplasmicMembrane)\ndefault = CytoplasmicMembrane')
parser.add_argument('--output', '-o', type=str, default="", required=True, help='Output directory')






args = parser.parse_args()
file_name =  args.sequence
file_title = file_name.split(".")
if file_title[1] != 'fasta':
	sys.exit('Please provide a FASTA file')

file_input = file_title[0]

bac_type = args.type
if args.output == "":
	sys.exit('please provide a output directory')


os.system('mkdir temp > /dev/null')
# In[4]:
os.system("rm -r temp/psortb")

###Removal of Shorter Sequences
print("Removing shorter sequences")
rmv_short_seq = "seqtk seq -L {} {}.fasta > temp/{}_large_seq.fasta".format(args.length,file_input,file_input)
os.system(rmv_short_seq)





###Removal of Redundant Sequences
cd_hit = 'cd-hit -i temp/{}_large_seq.fasta -o temp/{}_non_redundant.fasta -c {} -n 3 > /dev/null'.format(file_input,file_input,args.cutoff)
os.system(cd_hit)

# In[6]:


###blast of processed sequences with the DEG database 
print("blast with DEG")
blast_deg_database = "blastp -query temp/{}_non_redundant.fasta -db database/DEG -out temp/query_DEG.csv -outfmt '6 qseqid qacc qlen sseqid sacc slen qstart qend sstart send evalue bitscore length qcovs pident'  -evalue 1e-100 -max_target_seqs 1 > /dev/null".format(file_input)
os.system(blast_deg_database)
#df_path = pd.read_csv("temp/query_pathogen.csv", sep = "\t", names = ['query' , 'subject', 'identity', "Alignment lenth", "mismatches", "gap opens", "qstart", "qend", "s start", "s end", "evalue", "bitscore"])
df_deg = pd.read_csv("temp/query_DEG.csv", sep = "\t", names = ["qseqid", "qlen", "sseqid", "sacc", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qcovs", "pident "] )
df_deg = df_deg.drop_duplicates(subset='qseqid', keep="last")


# In[7]:


df1_deg = df_deg
#df1_deg = df1_deg[df1_path.qcovs >= 95]
#df1_deg = df1_deg[df1_path.pident >= 85]
df1_deg = df1_deg['qseqid']

df1_deg.to_csv('temp/query_deg.list', index = False, header = False)
#print(df1_path)
seqtk1 = "seqtk subseq {}.fasta temp/query_deg.list > temp/deg_hit.fasta".format(file_input)
os.system(seqtk1)



###blast of processed sequences with the human database 
print("blast with human")
blast_human_database = "blastp -query temp/deg_hit.fasta -db database/human -out temp/query_human.csv -outfmt '6 qseqid qacc qlen sseqid sacc slen qstart qend sstart send evalue bitscore length qcovs pident'  -evalue 1e-4 -max_target_seqs 1 > /dev/null ".format(file_input)
os.system(blast_human_database)
#df_path = pd.read_csv("temp/query_pathogen.csv", sep = "\t", names = ['query' , 'subject', 'identity', "Alignment lenth", "mismatches", "gap opens", "qstart", "qend", "s start", "s end", "evalue", "bitscore"])
df_human = pd.read_csv("temp/query_human.csv", sep = "\t", names = ["qseqid", "qlen", "sseqid", "sacc", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qcovs", "pident "] )
df_human = df_human.drop_duplicates(subset='qseqid', keep="last")

df1_human = df_human
#df1_human = df1_human[df1_path.qcovs >= 95]
#df1_human = df1_human[df1_path.pident >= 85]
df1_human = df1_human['qseqid']

df1_human.to_csv('temp/query_human.list', index = False, header = False)
#print(df1_path)
seqtk2 = "seqtk subseq {}.fasta temp/query_human.list > temp/human_hit.fasta".format(file_input)
os.system(seqtk2)


human_hit_ids = set([rec.id for rec in SeqIO.parse('temp/human_hit.fasta', 'fasta')])
# get entries unique to pathogen_hit.fasta
non_human = (rec for rec in SeqIO.parse('temp/deg_hit.fasta', 'fasta') if rec.id not in human_hit_ids)
# write unique entries to new file: pathogen_unique.fna
with open('temp/non_human.fasta', 'w') as target:
    SeqIO.write(non_human, target, 'fasta')




###########make db




###blast of processed sequences with the kegg database 
print("blast with kegg")
blast_kegg_database = "blastp -query temp/non_human.fasta -db database/kegg -out temp/query_kegg.csv -outfmt '6 qseqid qacc qlen sseqid sacc slen qstart qend sstart send evalue bitscore length qcovs pident'  -evalue 1e-100 -max_target_seqs 1 > /dev/null".format(file_input)
os.system(blast_kegg_database)
#df_kegg = pd.read_csv("temp/query_kegg.csv", sep = "\t", names = ['query' , 'subject', 'identity', "Alignment lenth", "mismatches", "gap opens", "qstart", "qend", "s start", "s end", "evalue", "bitscore"])
df_kegg = pd.read_csv("temp/query_kegg.csv", sep = "\t", names = ["qseqid", "qlen", "sseqid", "sacc", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qcovs", "pident "] )
df_kegg = df_kegg.drop_duplicates(subset='qseqid', keep="last")

df1_kegg = df_kegg
#df1_kegg = df1_kegg[df1_kegg.qcovs >= 95]
#df1_kegg = df1_kegg[df1_kegg.pident >= 90]
df1_kegg = df1_kegg['qseqid']

df1_kegg.to_csv('temp/query_kegg.list', index = False, header = False)
#print(df1_kegg)
seqtk3 = "seqtk subseq {}.fasta temp/query_kegg.list > temp/kegg_hit.fasta".format(file_input)
os.system(seqtk3)

os.system("mkdir temp/psortb > /dev/null")

#make predictions for localizzation
psortb = "./psortb -{} -r temp/psortb -i temp/kegg_hit.fasta".format(bac_type)
os.system(psortb)
os.system("mv temp/psortb/*.txt temp/psortb/psortb_results.txt")


f = open("temp/psortb/psortb_results.txt", 'r')
rr = f.read()


outputdir = "mkdir {}".format(args.output)
os.system(outputdir)
protein = re.findall(r'SeqID:\s*(.*)\n', str(rr))
#print(protein)
df_loc = pd.DataFrame(data = protein)
df_loc.columns = ['Protein']
loc = []
location = re.findall(r'Final Prediction:\n(.*)\n',str(rr))
for i in location:
    k = re.sub(r'[0-9]+', '', str(i))
    k = k.replace('.', '')
    k = k.replace(' ', '')
    loc.append(k)
    

#print(loc)
df_loc["Localization"] = pd.Series(loc)
#print(loc)
#print(df_loc)
df_loc.to_csv('temp/localization.csv')

os.system("rm -r temp/psortb")

df_loc = df_loc[df_loc.Localization == args.localization]
df_loc = df_loc['Protein']
df_loc.to_csv('temp/selected_proteins.list', index = False, header = False)

denoise = open("temp/selected_proteins.list", 'r')
den = denoise.read()
den = den.replace('''"''', "")
with open("temp/selected_proteins.list", "w") as wrt:
    for idss in den:
        wrt.write(idss)


#print(df_loc)
seqtk4 = "seqtk subseq {}.fasta temp/selected_proteins.list > {}/selected_proteins.fasta".format(file_input,args.output)
protein_targets = "cp temp/selected_proteins.list {}/selected_proteins.txt".format(args.output)
os.system(seqtk4)
os.system(protein_targets)





# Go pipeline, go!!!
#if __name__ == "__main__":
#    main()


