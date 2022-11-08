# HIP_figures

Cada carpeta corresponde a un conjunto de genomas:
- Genomas completos
- Genomas completos y genomas solo con cromosoma
- Genomas usados en el estudio de elhai (2015)
- Genomas usados en el estudio  de PJ Cabello-Yeves (2022). Incluyen53 nuevos genomas de picocyanobacterias.

Cada conjunto de genomas tiene 4 conjuntos de analisis de conteos de palindromos que corresponden con 4 valores de significancia que se tomaron como umbral de acuerdo al qvalue del FDR: 1e-32, 1e-64, 1e-128, 1e-256.

Para cada conjunto de analisis hay una grafica de puntos que muestra aquellos PALINDROMOS SIGNIFICATIVOS para el el umbral de corte indicado (1e-32, 1e-64, 1e-128 ó 1e-256).

Para cada conjunto de Genomas se hicieron 6 analisis de K-meros:
- Octanucleotidos (Octnaucs)
- Hexanucleotidos (Hexanucs)
- Pentanucleotidos (Pentanucs)
- Tetranucleotidos (Tetranucs)
- Hexanucleotidos sin contexto en HIP1 (HexanucsNoHIP)
- Tetranucleotidos sin contexto en HIP1 (TetranucsNoHIP)

Para cada analisis de K-meros se anotaron 2 filogenias de acuerdo a la abundancia:
- Observado/Esperado (OE)
- Frecuencia Observada cada 1000 nt (FrecObs)

Cada arbol contiene:
- Un Heatmap que muestra las abundancias significativas de todos los palindromos para cada especie y
- Un Barplot que muestra el palindromo mas abundante de entre todos los significativos para esa especie.

Tambien se incluye el Script de R usado para anotar las filogenias.

