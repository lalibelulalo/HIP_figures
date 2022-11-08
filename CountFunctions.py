from Bio import SeqIO
import re
import pandas as pd
import numpy as np

def Markov(kmer,genome_length,order,NUC,DINUC,TRINUC,TETRANUC):
    order = int(order)
    if (order == 0):
        K1NUCLEOTIDES = NUC
    elif (order == 1):
        K1NUCLEOTIDES = DINUC
        K2NUCLEOTIDES = NUC
    elif (order == 2):
        K1NUCLEOTIDES = TRINUC
        K2NUCLEOTIDES = DINUC
    elif (order == 3):
        K1NUCLEOTIDES = TETRANUC
        K2NUCLEOTIDES = TRINUC    
    if (order != 0):
        ## MODELO DE ORDEN K:
        # Producto de frecuencias de [K]nucleotidos del palindromo en el genoma (NUMERADOR)
        numerador_frec_K1nuc = 1
        for ii, ch in enumerate(kmer):
            K1nuc = kmer[ii:ii + (order+1)]
            if not (len(K1nuc) == (order+1)):
                break
            numerador_frec_K1nuc = K1NUCLEOTIDES[K1nuc] * (1 / (genome_length - order)) * numerador_frec_K1nuc
        # Producto de frecuencias de [K-1]nucleotidos del palindromo en el genoma (DENOMINADOR)
        contador_extremos_K2nuc, denominador_frec_K2nuc = 1, 1
        for ii, ch in enumerate(kmer):
            K2nuc = kmer[ii + 1:ii + (order+1)]
            if (contador_extremos_K2nuc <= (len(kmer)-(order-1))-2):
                contador_extremos_K2nuc += 1
                if not (len(K2nuc) == order):
                    break
                denominador_frec_K2nuc = K2NUCLEOTIDES[K2nuc] * (1 / (genome_length - (order-1))) * denominador_frec_K2nuc

        # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN K)
        if not (denominador_frec_K2nuc != 0):
            markov = 0
        else:
            markov = (numerador_frec_K1nuc / denominador_frec_K2nuc) * (genome_length - len(kmer)-1)
    else:
        numerador_frec_K1nuc = 1
        ## MODELO DE ORDEN 0:
        # Producto de frecuencias de nucleotidos del palindromo en el genoma (DENOMINADOR)
        frec_K1nuc = 1
        for ii, ch in enumerate(kmer):
            K1nuc = kmer[ii:ii + (order+1)]
            if not (len(nucleotide) == (order+1)):
                break
            frec_K1nuc = K1NUCLEOTIDES[K1nuc] * (1 / genome_length) * frec_K1nuc
            # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 0)
        markov = (frec_K1nuc) * (genome_length - len(kmer)-1)
    return(markov)

def Markov46(kmer,genome_length,order,NUC,DINUC,TRINUC,TETRANUC,HONUC):
    order = int(order)
    if kmer == 'GATC':
        pal='CGATCG'
        try:
            EX = HONUC[pal]
        except KeyError:
            EX = 0
    if kmer == 'CGATCG':
        pal='GCGATCGC'
        try:
            EX = HONUC[pal]
        except KeyError:
            EX = 0
        
    if (order == 0):
        K1NUCLEOTIDES = NUC
    elif (order == 1):
        K1NUCLEOTIDES = DINUC
        K2NUCLEOTIDES = NUC
    elif (order == 2):
        K1NUCLEOTIDES = TRINUC
        K2NUCLEOTIDES = DINUC
    elif (order == 3):
        K1NUCLEOTIDES = TETRANUC
        K2NUCLEOTIDES = TRINUC    
    if (order != 0):
        ## MODELO DE ORDEN K:
        # Producto de frecuencias de [K]nucleotidos del palindromo en el genoma (NUMERADOR)
        numerador_frec_K1nuc = 1
        for ii, ch in enumerate(kmer):
            K1nuc = kmer[ii:ii + (order+1)]
            if not (len(K1nuc) == (order+1)):
                break
            if K1nuc == 'CG' or K1nuc == 'GC':
                numerador_frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX*2)) * (1 / (genome_length - order)) * numerador_frec_K1nuc
            else:
                numerador_frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX)) * (1 / (genome_length - order)) * numerador_frec_K1nuc
                

        # Producto de frecuencias de [K-1]nucleotidos del palindromo en el genoma (DENOMINADOR)
        contador_extremos_K2nuc, denominador_frec_K2nuc = 1, 1
        for ii, ch in enumerate(kmer):
            K2nuc = kmer[ii + 1:ii + (order+1)]
            if (contador_extremos_K2nuc <= (len(kmer)-(order-1))-2):
                contador_extremos_K2nuc += 1
                if not (len(K2nuc) == order):
                    break
                if K1nuc == 'CG' or K1nuc == 'GC':
                    denominador_frec_K2nuc = (K2NUCLEOTIDES[K2nuc]-(EX*2)) * (1 / (genome_length - (order-1))) * denominador_frec_K2nuc
                else:
                    denominador_frec_K2nuc = (K2NUCLEOTIDES[K2nuc]-(EX)) * (1 / (genome_length - (order-1))) * denominador_frec_K2nuc

        # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN K)
        if not (denominador_frec_K2nuc != 0):
            markov = 0
        else:
            markov = (numerador_frec_K1nuc / denominador_frec_K2nuc) * (genome_length - len(kmer)-1)
    else:
        numerador_frec_K1nuc = 1
        ## MODELO DE ORDEN 0:
        # Producto de frecuencias de nucleotidos del palindromo en el genoma (DENOMINADOR)
        frec_K1nuc = 1
        for ii, ch in enumerate(kmer):
            K1nuc = kmer[ii:ii + (order+1)]
            if not (len(nucleotide) == (order+1)):
                break
            if pal == 'CGATCG':
                if K1nuc == 'C' or K1nuc == 'G':
                    frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX*2)) * (1 / genome_length) * frec_K1nuc
                else:
                    frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX)) * (1 / genome_length) * frec_K1nuc
            if pal == 'GCGATCGC':
                if K1nuc == 'C' or K1nuc == 'G':
                    frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX*3)) * (1 / genome_length) * frec_K1nuc
                else:
                    frec_K1nuc = (K1NUCLEOTIDES[K1nuc]-(EX)) * (1 / genome_length) * frec_K1nuc
                
                
            # DIVISION ENTRE FRACCIONES (MODELO DE MARKOV ORDEN 0)
        markov = (frec_K1nuc) * (genome_length - len(kmer)-1)
    return(markov)

def Gbk2Fna(GenomeDir,GenomeFile):
    gbk_file = str("".join ([GenomeDir,GenomeFile]))
    fna_filename = re.sub('_genomic.gbff', '.fna', GenomeFile)
    input_handle  = open(gbk_file, "r")
    output_handle = open(fna_filename, "w")
    for seq_record in SeqIO.parse(input_handle, "genbank") :
        output_handle.write(">%s %s\n%s\n" % (
            seq_record.id,
            seq_record.description,
            seq_record.seq))
    output_handle.close()
    input_handle.close()
    return(fna_filename)

def NucsGenLen(genome):
    # BUSCAMOS NUCLEOTIDOS EN TODO EL GENOMA Y OBTENEMOS EL TAMAÑO DEL GENOMA
    nuc_genome, nuc_list = [genome[ii:ii+1] for ii, ch in enumerate(genome)], ('A','T','C','G','N')
    NUCLEOTIDES, genome_length = {nuc:nuc_genome.count(nuc) for nuc in nuc_list}, len(nuc_genome)
    return(NUCLEOTIDES,genome_length)

def CountNucsAndKmers(genome,kmer):
    # k-meros del genoma
    if kmer >= 4:
        dinuc_genome, trinuc_genome, tetranuc_genome = [genome[ii:ii+2] for ii, ch in enumerate(genome)], [genome[ii:ii+3] for ii, ch in enumerate(genome)], [genome[ii:ii+4] for ii, ch in enumerate(genome)]
        dinuc_genome_cnts, trinuc_genome_cnts, tetranuc_genome_cnts = pd.Series(dinuc_genome).value_counts().to_dict(), pd.Series(trinuc_genome).value_counts().to_dict(), pd.Series(tetranuc_genome).value_counts().to_dict()
        NucsKmers = (dinuc_genome_cnts, trinuc_genome_cnts, tetranuc_genome_cnts)
        
        if kmer >= 5:
            pentanuc_genome = [genome[ii:ii+5] for ii, ch in enumerate(genome)]
            pentanuc_genome_cnts = pd.Series(pentanuc_genome).value_counts().to_dict()
            NucsKmers = (dinuc_genome_cnts, trinuc_genome_cnts, tetranuc_genome_cnts, pentanuc_genome_cnts)
            
            if kmer >= 6 or kmer == 41:
                hexanuc_genome = [genome[ii:ii+6] for ii, ch in enumerate(genome)]
                hexanuc_genome_cnts = pd.Series(hexanuc_genome).value_counts().to_dict()
                NucsKmers = (dinuc_genome_cnts, trinuc_genome_cnts, tetranuc_genome_cnts, pentanuc_genome_cnts, hexanuc_genome_cnts)
                
                if kmer == 8 or kmer == 82 or kmer == 61:
                    octanuc_genome = [genome[ii:ii+8] for ii, ch in enumerate(genome)]
                    octanuc_genome_cnts = pd.Series(octanuc_genome).value_counts().to_dict()
                    NucsKmers = (dinuc_genome_cnts, trinuc_genome_cnts, tetranuc_genome_cnts, pentanuc_genome_cnts, hexanuc_genome_cnts, octanuc_genome_cnts)
                    
        return(NucsKmers)

def PalKmerCounts(pal,kmer,NucsAndKmers):
    if kmer >= 4:
        ## Dinucleotidos del palindromo en cuestion
        dinuc_genome_cnts = NucsAndKmers[0]
        dinuc_pal = [pal[ii:ii+2] for ii, ch in enumerate(pal) if len(pal[ii:ii+2])==2 ]
        DINUCLEOTIDES = {x:(dinuc_genome_cnts[x] if x in dinuc_genome_cnts else 0) for x in dinuc_pal}

        ## Trinucleotidos del palindromo en cuestion
        trinuc_genome_cnts = NucsAndKmers[1]
        trinuc_pal = [pal[ii:ii+3] for ii, ch in enumerate(pal) if len(pal[ii:ii+3])==3 ]
        TRINUCLEOTIDES = {x:(trinuc_genome_cnts[x] if x in trinuc_genome_cnts else 0) for x in trinuc_pal}

        ## Tetranucleotidos del palindromo en cuestion
        tetranuc_genome_cnts = NucsAndKmers[2]
        tetranuc_pal = [pal[ii:ii+4] for ii, ch in enumerate(pal) if len(pal[ii:ii+4])==4 ]
        TETRANUCLEOTIDES = {x:(tetranuc_genome_cnts[x] if x in tetranuc_genome_cnts else 0) for x in tetranuc_pal}
        KNUCLEOTIDES = (DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES)
        
        if kmer >= 5:
            pentanuc_genome_cnts = NucsAndKmers[3]
            ## Pentanucleotidos del palindromo en cuestion
            pentanuc_pal = [pal[ii:ii+5] for ii, ch in enumerate(pal) if len(pal[ii:ii+5])==5 ]
            PENTANUCLEOTIDES = {x:(pentanuc_genome_cnts[x] if x in pentanuc_genome_cnts else 0) for x in pentanuc_pal}
            KNUCLEOTIDES = (DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES)
        
            if kmer >= 6 or kmer == 41:
                hexanuc_genome_cnts = NucsAndKmers[4]
                ## Hexanucleotidos del palindromo en cuestion
                hexanuc_pal = [pal[ii:ii+6] for ii, ch in enumerate(pal) if len(pal[ii:ii+6])==6 ]
                HEXANUCLEOTIDES = {x:(hexanuc_genome_cnts[x] if x in hexanuc_genome_cnts else 0) for x in hexanuc_pal}
                KNUCLEOTIDES = (DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES,HEXANUCLEOTIDES)
                
                if kmer == 8 or kmer == 82 or kmer == 61:
                    octanuc_genome_cnts = NucsAndKmers[5]
                    ## Octanucleotidos del palindromo en cuestion
                    octanuc_pal = [pal[ii:ii+8] for ii, ch in enumerate(pal) if len(pal[ii:ii+8])==8 ]
                    OCTANUCLEOTIDES = {x:(octanuc_genome_cnts[x] if x in octanuc_genome_cnts else 0) for x in octanuc_pal}
                    KNUCLEOTIDES = (DINUCLEOTIDES,TRINUCLEOTIDES,TETRANUCLEOTIDES,PENTANUCLEOTIDES,HEXANUCLEOTIDES,OCTANUCLEOTIDES)
    return(KNUCLEOTIDES)
    
def TetranucGen(Nucs):
    pals = [Nucs[i]+Nucs[j]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    return(pals)

def TetranucNoHIPGen(Nucs):
    pals = [Nucs[i]+"GATC"+Nucs[j] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    pals.remove('CGATCG')
    return(pals)

def PentanucGen(Nucs):
    pals = [Nucs[i]+Nucs[j]+Nucs[k]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs))) for k in range(1,(len(Nucs)))]
    return(pals)
    
def HexanucGen(Nucs):
    pals = [Nucs[i]+Nucs[j]+Nucs[k]+Nucs[-k]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs))) for k in range(1,(len(Nucs)))]
    return(pals)

def HexanucNoHIPGen(Nucs):
    pals = ["".join ([Nucs[i],"CGATCG",Nucs[j]]) for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    pals.remove('GCGATCGC')
    return(pals)

def OctanucGen(Nucs):
    pals = [Nucs[i]+Nucs[j]+Nucs[k]+Nucs[l]+Nucs[-l]+Nucs[-k]+Nucs[-j]+Nucs[-i] for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs))) for k in range(1,(len(Nucs))) for l in range(1,(len(Nucs)))]
    return(pals)

def OctanucHIPGen(Nucs):
    hip_like_1 = ["".join ([Nucs[i],"CGATCG",Nucs[j]]) for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    hip_like_2 = ["G"+Nucs[i]+"GATC"+Nucs[j]+"C" for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    hip_like_3 = ["GC"+Nucs[i]+"AT"+Nucs[j]+"GC" for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    hip_like_4 = ["GCG"+Nucs[i]+Nucs[j]+"CGC" for i in range(1,(len(Nucs))) for j in range(1,(len(Nucs)))]
    hip_like_all = [hip_like_1,hip_like_2,hip_like_3,hip_like_4]
    hip_like_all = [j for i in hip_like_all for j in i]
    hip_like = np.unique(hip_like_all)
    pals = hip_like
    return(pals)
