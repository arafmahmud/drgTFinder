
import io
import os
import re
from Bio import SeqIO
from Bio.KEGG import REST
import pandas as pd
import argparse




parser = argparse.ArgumentParser(prog="make_db.py", formatter_class=argparse.RawTextHelpFormatter, description="""
make database for the drgTFinder  \n""",epilog="""
    Examples:
    python3 make_db.py -b eco -h GRCH38.fasta -d DEG.fasta\n""")
parser.add_argument('--bacteria', '-b', type=str, default="", required=True, help='Three or four letter code for the query organism found in https://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=show.org')
parser.add_argument('--host', '-hs', type=str, default="", required=True, help='Host whole proteome as fasta file found in https://www.uniprot.org/proteomes/UP000005640')
parser.add_argument('--DEG', '-d', type=str, default="", required=True, help='DEG database fasta file found in http://tubic.tju.edu.cn/deg/download.php')

args = parser.parse_args()


os.system("mkdir database")
hm = REST.kegg_list("pathway", "hsa").read()
hm = hm.replace(" - Homo sapiens (human)", "")
with open("database/human_pathway.csv", "w") as f1:
    for j in hm:
        f1.write(j)
df2 = pd.read_csv("database/human_pathway.csv", sep = "\t", names = ["Accession", "Pathway"])

gene_list = REST.kegg_list("bmf").read() #.format()
with open("database/gene_list.csv", "w") as f:
    f.write(gene_list)
    
df_gene_list = pd.read_csv("database/gene_list.csv", sep = "\t", names = ["Accession", "Name"])
#print(df)
f2 = open("database/unique_pathway_proteins.fasta", "w") 

for ids in df_gene_list["Accession"]:
    #print(ids)

    get_info = REST.kegg_get(str(ids)).read()
    
    pathway = re.findall(r'PATHWAY\s*(.*)MODULE', str(get_info), re.DOTALL)
    if pathway == []:
        pass
    else:
        print(pathway)
#    
        with open("database/pathway.csv", "w") as f:
            for i in pathway:
                s = re.sub("[0-9](\s+)[A-Z]", "\t", i)
                f.write(s)
        df_pathway = pd.read_csv("database/pathway.csv", sep = "\t", names = ["Accession", "Pathway"])
        #print(df)
        k = df_pathway["Pathway"].isin(df2["Pathway"])
        k = str(k)
        if "True" in k:
            pass
        else:
            print("Unique!!")
            r = REST.kegg_get(ids, "aaseq").read()
            f2.write(r)
            print(r)
            
#
#
make_kegg_database = "makeblastdb -in database/unique_pathway_proteins.fasta -out database/kegg  -dbtype prot"
os.system(make_kegg_database)
make_human_database = "makeblastdb -in {} -out database/human  -dbtype prot".format(args.host)
os.system(make_human_database)
make_deg_database = "makeblastdb -in {} -out database/DEG  -dbtype prot".format(args.DEG)
os.system(make_deg_database)
#print(df2)

#

##        
#
#'''
#    '''
#    
