# drgTFinder
drgTFinder in a automated tool for finding Drug Targets from bacterial whole proteome sequence.

# Prerequisites:

Create a new conda environment and install dependencies:
```
conda create --name drgTFinder 
conda activate drgTFinder
conda install -c bioconda -y blast
conda install -c bioconda -y seqtk
pip3 install biopython
pip3 install pandas

```

# Make Database
```
make_db.py [-h] --bacteria BACTERIA --host HOST --DEG DEG
                  [--database DATABASE] [--threads THREADS]

make database for the drgTFinder  

optional arguments:
  -h, --help            show this help message and exit
  --bacteria BACTERIA, -b BACTERIA
                        Three or four letter code for the query organism found in https://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=show.org
  --host HOST, -hs HOST
                        Host whole proteome as fasta file found in https://www.uniprot.org/proteomes/UP000005640
  --DEG DEG, -d DEG     DEG database fasta file found in http://tubic.tju.edu.cn/deg/download.php
  --database DATABASE, -db DATABASE
                        Database directory
  --threads THREADS, -t THREADS
                        Number of threads

    Examples:
    python3 make_db.py -b cvi -hs GRCH38.fasta -d DEG.fasta -t 8 -db database 
```
# Use drgTFinder:
```
usage: drgTFinder.py [-h] --sequence SEQUENCE [--length LENGTH]
                     [--cutoff CUTOFF] [--output OUTPUT] [--database DATABASE]
                     [--threads THREADS]

drgTFinder is an automated pipeline for detecting Novel drug targets using proteome sequence from 

optional arguments:
  -h, --help            show this help message and exit
  --sequence SEQUENCE, -s SEQUENCE
                        Whole proteome sequence file in FASTA format.
  --length LENGTH, -l LENGTH
                        Sequence Lenth cutoff (default = 100)
  --cutoff CUTOFF, -c CUTOFF
                        CD-Hit cutoff value (default = 0.8)
  --output OUTPUT, -o OUTPUT
                        Output directory
  --database DATABASE, -db DATABASE
                        Database directory
  --threads THREADS, -t THREADS
                        Number of threads

    Examples:
    python3 drgTFinder.py -db database -s test.fasta -o output -t 4

```
