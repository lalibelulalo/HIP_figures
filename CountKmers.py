from Bio import SeqIO
from CountFunctions import Gbk2Fna
from CountFunctions import NucsGenLen
from CountFunctions import CountNucsAndKmers
from CountFunctions import PalKmerCounts
from CountFunctions import Markov
from CountFunctions import Markov46
from CountFunctions import TetranucGen
from CountFunctions import TetranucNoHIPGen
from CountFunctions import PentanucGen
from CountFunctions import HexanucGen
from CountFunctions import HexanucNoHIPGen
from CountFunctions import OctanucGen
from CountFunctions import OctanucHIPGen
import os
import glob
import time
import re
    
def CountKmers(GENOMES_PATH,KMERS,ORDER):
    ## Obtenemos tiempo inicial
    st = time.time()
    st2 = time.time()
    DateTime = time.strftime("%Y,%m,%d,%H,%M,%S")
    t = DateTime.split(',')
    numbers = [ int(x) for x in t ]  
    current_date = "-".join([str(numbers[0]),str(numbers[1]),str(numbers[2])])
    current_time = "".join ([str(numbers[3]),'hrs',str(numbers[4]),'mins'])
    date = str("_".join ([current_date,current_time]))
    ##-----------------------------------------------
    ## CARGAMOS PALINDROMOS
    Nucs = ["N","A","G","C","T"]
    ## Octameros
    if KMERS == 4:
        pals = TetranucGen(Nucs)
        kmer = 'Tetranuc'
    elif KMERS == 41:
        pals = TetranucNoHIPGen(Nucs)
        kmer = 'TetranucNoHIP'
    elif KMERS == 5:
        pals = PentanucGen(Nucs)
        kmer = 'Pentanuc'
    elif KMERS == 6:
        pals = HexanucGen(Nucs)
        kmer = 'Hexanuc'
    elif KMERS == 61:
        pals = HexanucNoHIPGen(Nucs)
        kmer = 'HexanucNoHIP'
    elif KMERS == 8:
        pals = OctanucGen(Nucs)
        kmer = 'Octanuc'
    elif KMERS == 82:
        pals = OctanucHIPGen(Nucs)
        kmer = 'OctanucHIP'
    else:
        print("Opcion incorrecta")
    ##-----------------------------------------------    
    ## Nombre del modelo
    model = str("".join(['M',str(ORDER)]))
    ## Nombre del archivo de salida
    if GENOMES_PATH.endswith("/"):
        genomepath = re.sub('/', '', GENOMES_PATH)
    else:
        genomepath = GENOMES_PATH
    output_file = str("_".join (['Markov_count',genomepath,date,kmer,model,'.txt']))

    ## Imprimimos el nombre del archivo de salida
    print ("The output file is: {}.\n".format(output_file))

    ## Creamos el archivo de salida
    header = str("".join(['spp\tID\tpalindrome\tobs\tmarkov',str(ORDER),'\tgenomesize\tA\tTh\tC\tG\tN\n']))
    output = open (output_file, 'w')
    output.write(header)

    ## Entramos al path proporcionado que contiene los genomas
    if GENOMES_PATH.endswith("/"):
        GenomeDir = str(GENOMES_PATH)
    else:
        GenomeDir = str("".join ([GENOMES_PATH,'/']))
    ## Hacemos una lista con los genomas
    genomes = [x for x in os.listdir(GenomeDir) if x.endswith(".gbff") or x.endswith(".gbk")]

    contador= 0
    for GenomeFile in genomes:
        KNUCSUM = 0
        contador += 1
        print ("Archivo {} de {}: {}".format(contador,len(genomes), GenomeFile))
        ## Obtenemos nombre del genoma en turno
        GenomeFileName = Gbk2Fna(GenomeDir,GenomeFile)
        GenomeFile = GenomeFileName #Gbk2Fna(GenomeDir,GenomeFile)     
        ## Obtenemos el ID
        try:
            ID = re.search(r'GCF_\d+.\d+', GenomeFile)
            if ID is None:
                ID = re.search(r'GCA_\d+.\d+', GenomeFile)
                if ID is None:
                    ID = re.search(r'\w\w_\d+.\d+', GenomeFile)
            ID = ID.group()

        except AttributeError:
            ID = 'GCX_000000000.0'

        NameEnd = str("".join(["_",str(ID),".fna"]))#.fna
        Spp = GenomeFileName.replace(NameEnd, "")
        ## Extraemos las secuencias sin encabezados
        fasta_sequences = SeqIO.parse(open(GenomeFile),'fasta')
        genome = [str(fasta.seq) for fasta in fasta_sequences]
        ## Creamos un string con todo el genoma
        genome = ''.join(genome)
        genome = genome.upper()

        # BUSCAMOS NUCLEOTIDOS EN TODO EL GENOMA Y OBTENEMOS EL TAMAÑO DEL GENOMA
        NUCLEOTIDES, genome_length = NucsGenLen(genome)

        # k-meros del genoma
        NucsKmers = CountNucsAndKmers(genome,KMERS)
            
        for pal in pals:
            pal=pal.upper()
            if KMERS == 4:
                DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES = PalKmerCounts(pal,KMERS,NucsKmers)
                KNUCS = TETRANUCLEOTIDES
            elif KMERS == 5:
                DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES = PalKmerCounts(pal,KMERS,NucsKmers)
                KNUCS = PENTANUCLEOTIDES
            elif KMERS == 6 or KMERS == 41:
                DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES,HEXANUCLEOTIDES = PalKmerCounts(pal,KMERS,NucsKmers)
                KNUCS = HEXANUCLEOTIDES
            elif KMERS == 8 or KMERS == 82 or KMERS == 61:
                DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES,HEXANUCLEOTIDES,OCTANUCLEOTIDES = PalKmerCounts(pal,KMERS,NucsKmers)
                KNUCS = OCTANUCLEOTIDES
            else:
                print("Opcion incorrecta")
            if KMERS == 4 or KMERS == 5 or KMERS == 6 or KMERS == 8 or KMERS == 82:
                markov = Markov(pal,genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES)
                output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Spp,ID,pal,KNUCS[pal],markov,genome_length,NUCLEOTIDES['A'],NUCLEOTIDES['T'],NUCLEOTIDES['C'],NUCLEOTIDES['G'],NUCLEOTIDES['N']))
            elif KMERS == 41:
                #markov = Markov('GATC',genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES)
                KNUCSUM = KNUCSUM + KNUCS[pal]
            elif KMERS == 61:
                #markov = Markov('CGATCG',genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES)
                KNUCSUM = KNUCSUM + KNUCS[pal]
        if KMERS == 41:
            KNUCS_G = NucsKmers[4]
            markov = Markov46('GATC',genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,KNUCS_G)
            output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Spp,ID,'GATC',KNUCSUM,markov,genome_length,NUCLEOTIDES['A'],NUCLEOTIDES['T'],NUCLEOTIDES['C'],NUCLEOTIDES['G'],NUCLEOTIDES['N']))
        if KMERS == 61:
            KNUCS_G = NucsKmers[5]
            markov = Markov46('CGATCG',genome_length,ORDER,NUCLEOTIDES,DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,KNUCS_G)
            output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(Spp,ID,'CGATCG',KNUCSUM,markov,genome_length,NUCLEOTIDES['A'],NUCLEOTIDES['T'],NUCLEOTIDES['C'],NUCLEOTIDES['G'],NUCLEOTIDES['N']))
        
        #os.remove(fna_filename)
        os.remove(GenomeFile)
        # get the end time
        et = time.time()
        # get the execution time
        elapsed_time2 = et - st2
        st2 = time.time()
        print('Execution time: {} seconds.'.format(elapsed_time2))
        print ("------------------------------------")
    output.close()

    elapsed_time = et - st
    print ('*** Execution time: {} mins. ***'.format(elapsed_time/60))