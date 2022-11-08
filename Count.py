## use:
## python Conteo.py genomes_dir KMERS MODEL
## Where:
## KMERS: 
##	A. Octamers
##	B. HIP OCTAMERS
##	C. HEXAMERS
## MODEL: 0,1,2,3
## python Conteo.py all_pico_2022_fna B 3
import sys
from CountKmers import CountKmers

try:
    GenomePath = sys.argv[1]
    kmers = sys.argv[2]
    Order = sys.argv[3]
    CountKmers(str(GenomePath),int(kmers),Order)
except FileNotFoundError:
    print ("There directory path \"{}\"is incorrect".format(GenomePath))
