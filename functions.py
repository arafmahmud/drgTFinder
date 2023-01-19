import os
import pandas as pd
import re
from Bio import SeqIO
from Bio.KEGG import REST
import concurrent.futures


def remove_short_sequeces(threshold, input, output):
    cmd =  "seqtk seq -L {} {} > {}".format(threshold, input, output)
    os.system(cmd)

def remove_redundant_sequeces(input,output,threshold):
    cmd = 'cd-hit -i {} -o {} -c {} -n 3'.format(input,output,threshold)
    os.system(cmd)

def blast(input, output, database,threads):
    blast_run = "blastp -query {} -db {} -out {} -outfmt '6 qseqid qacc qlen sseqid sacc slen qstart qend sstart send evalue bitscore length qcovs pident'  -evalue 1e-100 -max_target_seqs 1 -num_threads {} > /dev/null".format(input,database, output, threads)
    os.system(blast_run)
def blast_host(input, output, database,threads):
    blast_run = "blastp -query {} -db {} -out {} -outfmt '6 qseqid qacc qlen sseqid sacc slen qstart qend sstart send evalue bitscore length qcovs pident'  -evalue 1e-4 -max_target_seqs 1 -num_threads {} > /dev/null".format(input,database, output, threads)
    os.system(blast_run) 


def generate_blast_hit_list(input, output):
    df = pd.read_csv(input, sep = "\t", names = ["qseqid", "qlen", "sseqid", "sacc", "slen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "qcovs", "pident "] )
    df = df.drop_duplicates(subset='qseqid', keep="last")
    df1 = df
    df1 = df1['qseqid']

    df1.to_csv(output, index = False, header = False)
    
    
def separate_hit_sequences(main_fasta, input, output):
    cmd = "seqtk subseq {} {} > {}".format(main_fasta, input, output)
    os.system(cmd)

def removal_of_sequences(input, separate, output):
    separate_ids = set([rec.id for rec in SeqIO.parse(separate, 'fasta')])
    # get entries unique to pathogen_hit.fasta
    input_fasta = (rec for rec in SeqIO.parse(input, 'fasta') if rec.id not in separate_ids)
    # write unique entries to new file: pathogen_unique.fna
    with open(output, 'w') as target:
        SeqIO.write(input_fasta, target, 'fasta')
    




def make_kegg_database(dbdir, organism,threads):
    os.system("mkdir {}".format(dbdir))
    hm = REST.kegg_list("pathway", "hsa").read()
    hm = hm.replace(" - Homo sapiens (human)", "")
    with open(dbdir + "/human_pathway.csv", "w") as f1:
        for j in hm:
            f1.write(j)
    df2 = pd.read_csv(dbdir + "/human_pathway.csv", sep = "\t", names = ["Accession", "Pathway"])

    #print(df2)
    gene_list = REST.kegg_list(organism).read() #.format()
    with open(dbdir + "/gene_list.tsv", "w") as f:
        f.write(gene_list)
        
    df_gene_list = pd.read_csv(dbdir + "/gene_list.tsv", sep = "\t", names = ["Accession","Type", "Location", "Name"])
    #print(df_gene_list)
    f2 = open(dbdir + "/unique_pathway_proteins.fasta", "w") 

    def multi(ids):
        print(ids)
        get_info = REST.kegg_get(str(ids)).read()
        print(get_info)
        
        pathway = re.findall(r'PATHWAY\s*(.*)MODULE', str(get_info), re.DOTALL)
        if pathway == []:
            pass
        else:
            print(pathway)
    #    
            with open(dbdir + "/pathway.csv", "w") as f:
                for i in pathway:
                    s = re.sub("[0-9](\s+)[A-Z]", "\t", i)
                    f.write(s)
            df_pathway = pd.read_csv(dbdir + "/pathway.csv", sep = "\t", names = ["Accession", "Pathway"])
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

    with concurrent.futures.ThreadPoolExecutor(max_workers = threads) as executor:
        run_preprocessor = [executor.submit(multi, accession) for accession in df_gene_list["Accession"]]

    
    
    # for ids in df_gene_list["Accession"]:
        
    #     print(ids)
    #     get_info = REST.kegg_get(str(ids)).read()
    #     print(get_info)
        
    #     pathway = re.findall(r'PATHWAY\s*(.*)MODULE', str(get_info), re.DOTALL)
    #     if pathway == []:
    #         pass
    #     else:
    #         print(pathway)
    # #    
    #         with open(dbdir + "/pathway.csv", "w") as f:
    #             for i in pathway:
    #                 s = re.sub("[0-9](\s+)[A-Z]", "\t", i)
    #                 f.write(s)
    #         df_pathway = pd.read_csv(dbdir + "/pathway.csv", sep = "\t", names = ["Accession", "Pathway"])
    #         #print(df)
    #         k = df_pathway["Pathway"].isin(df2["Pathway"])
    #         k = str(k)
    #         if "True" in k:
    #             pass
    #         else:
    #             print("Unique!!")
    #             r = REST.kegg_get(ids, "aaseq").read()
    #             f2.write(r)
    #             print(r)


def make_database(dbdir, kegg, host, deg):
    
    make_kegg_database = ''' cp {} {}/kegg.fasta
    makeblastdb -in {}/kegg.fasta -out {}/kegg  -dbtype prot'''.format(kegg, dbdir,dbdir,dbdir)

    make_host_database = ''' cp {} {}/host.fasta
    makeblastdb -in {}/host.fasta -out {}/host  -dbtype prot'''.format(host, dbdir,dbdir,dbdir)

    make_deg_database = ''' cp {} {}/deg.fasta
    makeblastdb -in {}/deg.fasta -out {}/deg  -dbtype prot'''.format(deg, dbdir,dbdir,dbdir)

    os.system(make_kegg_database)
    os.system(make_host_database)
    os.system(make_deg_database)
   
# def final_report(input):
#     separate_ids = set([rec.id for rec in SeqIO.parse(input, 'fasta')]) 
    
#     host_non_homologous = open("temp/host_homologous.list", "r")
#     deg_homologous = open( "temp/deg_homologous.list", "r")
#     kegg_non_homologous = open("temp/kegg_homologous.list", "r")

#     #print(f.read().splitlines())
#     #print(separate_ids)
#     #print(list(separate_ids))
#     for i in list(separate_ids):
#         #print(i)
#         if i not in host_non_homologous.read().splitlines():
#             print(i)
            
#         if i in deg_homologous.read().splitlines():
#             print(i)
#         if i in kegg_non_homologous.read().splitlines():
#             print(i)



#make_kegg_database("test_database", "cvi")
