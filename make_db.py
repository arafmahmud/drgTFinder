
import io
import os
import re
from Bio import SeqIO
from Bio.KEGG import REST
import pandas as pd
import argparse
from functions import make_database,make_kegg_database



parser = argparse.ArgumentParser(prog="make_db.py", formatter_class=argparse.RawTextHelpFormatter, description="""
make database for the drgTFinder  \n""",epilog="""
    Examples:
    python3 make_db.py -b cvi -hs GRCH38.fasta -d DEG.fasta -t 8 -db database \n""")
parser.add_argument('--bacteria', '-b', type=str, default="", required=True, help='Three or four letter code for the query organism found in https://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=show.org')
parser.add_argument('--host', '-hs', type=str, default="", required=True, help='Host whole proteome as fasta file found in https://www.uniprot.org/proteomes/UP000005640')
parser.add_argument('--DEG', '-d', type=str, default="", required=True, help='DEG database fasta file found in http://tubic.tju.edu.cn/deg/download.php')
parser.add_argument('--database', '-db', type=str, default="", required=False, help='Database directory')
parser.add_argument('--threads', '-t', type=int, default="4", required=False, help='Number of threads')




args = parser.parse_args()


os.system("mkdir {}".format(args.database))


make_kegg_database(args.database, args.bacteria, args.threads)
make_database(args.database, args.database+"/unique_pathway_proteins.fasta", args.host, args.DEG )





