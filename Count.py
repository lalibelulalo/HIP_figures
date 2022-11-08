## use:
## python Conteo.py genomes_dir KMERS MODEL
## Where:
## KMERS: 
##	8. Octamers
##	6. Hexamers
##  5. Pentamers
##  4. Tetramers
##  61. HexamersNoHIP
##	41. TetramersNoHIP
## MODEL: 0,1,2,3
## python Conteo.py all_pico_2022_fna 41 2
import sys
from CountKmers import CountKmers

try:
    GenomePath = sys.argv[1]
    kmers = sys.argv[2]
    Order = sys.argv[3]
    CountKmers(str(GenomePath),int(kmers),Order)
except FileNotFoundError:
    print ("There directory path \"{}\"is incorrect".format(GenomePath))
