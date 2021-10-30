# drgTFinder
Drug Target finder from bacterial whole proteome

Prerequisites:

Python libraries pandas,biopython are reuired installed via `pip3 install pandas biopython`

NCBI BLAST+ installed via `sudo apt-get install ncbi-blast+-legacy`

Psortb installed from https://hub.docker.com/r/brinkmanlab/psortb

# Make Database
`
usage: make_db.py [-h] --bacteria BACTERIA --host HOST --DEG DEG

make database for the drgTFinder  

optional arguments:
  -h, --help            show this help message and exit
  --bacteria BACTERIA, -b BACTERIA
                        Three or four letter code for the query organism found in https://www.genome.jp/kegg-bin/find_org_www?mode=abbr&obj=show.org
  --host HOST, -hs HOST
                        Host whole proteome as fasta file found in https://www.uniprot.org/proteomes/UP000005640
  --DEG DEG, -d DEG     DEG database fasta file found in http://tubic.tju.edu.cn/deg/download.php

    Examples:
    python3 make.py -b eco -h GRCH38.fasta -d DEG.fasta
`
# Use drgTFinder:
`
usage: drgTFinder.py [-h] --sequence SEQUENCE --type TYPE [--length LENGTH] [--cutoff CUTOFF] [--localization LOCALIZATION] --output OUTPUT

drgTFinder is an automated pipeline for detecting Novel drug targets using proteome sequence from 

optional arguments:
  -h, --help            show this help message and exit
  --sequence SEQUENCE, -s SEQUENCE
                        Whole proteome sequence file in FASTA format.
  --type TYPE, -t TYPE  Bacteria type(p or n) 
  --length LENGTH, -l LENGTH
                        Sequence Lenth cutoff (default = 100)
  --cutoff CUTOFF, -c CUTOFF
                        CD-Hit cutoff value (default = 0.8)
  --localization LOCALIZATION, -loc LOCALIZATION
                        Subcellular localization (Cytoplasmic, Extracellular, Cellwall, CytoplasmicMembrane)
                        default = CytoplasmicMembrane
  --output OUTPUT, -o OUTPUT
                        Output directory

    Examples:
    python3 drgTFinder.py -s whole_proteome.fasta -o outputdir -t n

`
