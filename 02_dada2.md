R Notebook
================
Pierre Goasdoue

  - [Dada2](#dada2)
      - [Importation des librairies et préparation des
        données](#importation-des-librairies-et-préparation-des-données)
      - [Inspection des profils de qualité des
        reads](#inspection-des-profils-de-qualité-des-reads)
      - [Filtration et tri des
        donnnées](#filtration-et-tri-des-donnnées)
      - [Calcul des erreurs et visualisation de ces
        derniers](#calcul-des-erreurs-et-visualisation-de-ces-derniers)
      - [Application de Dada2 aux
        données](#application-de-dada2-aux-données)
      - [Alignement des séquences Forward et
        Reverse](#alignement-des-séquences-forward-et-reverse)
      - [Table d’observation](#table-dobservation)
      - [Elimination des séquences
        chimériques](#elimination-des-séquences-chimériques)
      - [Résumé des opérations effectuées
        précédement](#résumé-des-opérations-effectuées-précédement)
      - [Téléchargement des données
        Silva](#téléchargement-des-données-silva)
  - [Phyloseq](#phyloseq)
      - [Importation des librairies et préparation des
        données](#importation-des-librairies-et-préparation-des-données-1)
      - [Visualisation de l’alpha
        diversité](#visualisation-de-lalpha-diversité)
      - [Ordination](#ordination)
      - [Représention des échantillons en “histogramme” (bar
        plot)](#représention-des-échantillons-en-histogramme-bar-plot)

\#1/ quelles sont les influences relative de la profondeur et de la
saison sur la structure des communautes planctoniques de la rade de
Brest \#2/ Quels sont les biomarkeurs de saison (hivers et ete) ?

# Dada2

``` r
library(Rcpp)
library(dada2)
```

    ## Warning: multiple methods tables found for 'which'

## Importation des librairies et préparation des données

``` r
path <- "~/CC2_Dada2AndPhyloseq/St_Stratif_11mars15"
list.files(path)
```

    ##  [1] "filtered"                            "Station5_Fond1_10sept14_R1.fastq"   
    ##  [3] "Station5_Fond1_10sept14_R2.fastq"    "Station5_Fond1_11mars15_R1.fastq"   
    ##  [5] "Station5_Fond1_11mars15_R2.fastq"    "Station5_Fond2_10sept14_R1.fastq"   
    ##  [7] "Station5_Fond2_10sept14_R2.fastq"    "Station5_Fond2_11mars15_R1.fastq"   
    ##  [9] "Station5_Fond2_11mars15_R2.fastq"    "Station5_Fond3_10sept14_R1.fastq"   
    ## [11] "Station5_Fond3_10sept14_R2.fastq"    "Station5_Median1_10sept14_R1.fastq" 
    ## [13] "Station5_Median1_10sept14_R2.fastq"  "Station5_Median2_10sept14_R1.fastq" 
    ## [15] "Station5_Median2_10sept14_R2.fastq"  "Station5_Surface1_10sept14_R1.fastq"
    ## [17] "Station5_Surface1_10sept14_R2.fastq" "Station5_Surface1_11mars15_R1.fastq"
    ## [19] "Station5_Surface1_11mars15_R2.fastq" "Station5_Surface2_10sept14_R1.fastq"
    ## [21] "Station5_Surface2_10sept14_R2.fastq" "Station5_Surface2_11mars15_R1.fastq"
    ## [23] "Station5_Surface2_11mars15_R2.fastq"

``` r
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "R"), `[`, 1)
```

## Inspection des profils de qualité des reads

``` r
plotQualityProfile(fnFs[1:2])
```

![](02_dada2_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

On observe une forte diminution du score de qualité vers 240 nucléotide.
La coupure se fera ici.

``` r
plotQualityProfile(fnRs[1:2])
```

![](02_dada2_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

On observe une forte diminution vers 200 nucléotides. Par rapport à
fnFs, le score de qualité est moins bien. La coupure se fera à 200.

## Filtration et tri des donnnées

``` r
filtFs <- file.path(path, "filtered", paste0(sample.names, "_Ffilt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_Rfilt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
```

``` r
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,200), trimLeft =21,
              maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=TRUE)
head(out)
```

    ##                                    reads.in reads.out
    ## Station5_Fond1_10sept14_R1.fastq     159971    145448
    ## Station5_Fond1_11mars15_R1.fastq     175993    160423
    ## Station5_Fond2_10sept14_R1.fastq     197039    177018
    ## Station5_Fond2_11mars15_R1.fastq      87585     79989
    ## Station5_Fond3_10sept14_R1.fastq     117140    106150
    ## Station5_Median1_10sept14_R1.fastq   116519    106745

## Calcul des erreurs et visualisation de ces derniers

``` r
errF <- learnErrors(filtFs, multithread=TRUE)
```

    ## 105752691 total bases in 482889 reads from 3 samples will be used for learning the error rates.

Estimation du taux d’erreur de filtFs

``` r
errR <- learnErrors(filtRs, multithread=TRUE)
```

    ## 100755162 total bases in 562878 reads from 4 samples will be used for learning the error rates.

Estimation du taux d’erreur de filtRs

``` r
plotErrors(errF, nominalQ=TRUE)
```

    ## Warning: Transformation introduced infinite values in continuous y-axis
    
    ## Warning: Transformation introduced infinite values in continuous y-axis

![](02_dada2_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Représentation des fréquences d’erreurs estimées. Points gris : taux
d’erreurs observées pour chacun des scores de qualité consensus. Ligne
noire : taux d’erreurs estimé après que l’algorithme ait réuni toutes
les informations liée aux taux d’erreurs estimés Ligne rouge : taux
d’erreurs attendu selon le Q-score

## Application de Dada2 aux données

``` r
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 37907 unique sequences.
    ## Sample 2 - 160423 reads in 35863 unique sequences.
    ## Sample 3 - 177018 reads in 47212 unique sequences.
    ## Sample 4 - 79989 reads in 20356 unique sequences.
    ## Sample 5 - 106150 reads in 30255 unique sequences.
    ## Sample 6 - 106745 reads in 28836 unique sequences.
    ## Sample 7 - 98823 reads in 25824 unique sequences.
    ## Sample 8 - 107427 reads in 26733 unique sequences.
    ## Sample 9 - 71082 reads in 17976 unique sequences.
    ## Sample 10 - 78645 reads in 20422 unique sequences.
    ## Sample 11 - 91534 reads in 24487 unique sequences.

``` r
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
```

    ## Sample 1 - 145448 reads in 45486 unique sequences.
    ## Sample 2 - 160423 reads in 41638 unique sequences.
    ## Sample 3 - 177018 reads in 55554 unique sequences.
    ## Sample 4 - 79989 reads in 23239 unique sequences.
    ## Sample 5 - 106150 reads in 34625 unique sequences.
    ## Sample 6 - 106745 reads in 31673 unique sequences.
    ## Sample 7 - 98823 reads in 29093 unique sequences.
    ## Sample 8 - 107427 reads in 28947 unique sequences.
    ## Sample 9 - 71082 reads in 21426 unique sequences.
    ## Sample 10 - 78645 reads in 22051 unique sequences.
    ## Sample 11 - 91534 reads in 28266 unique sequences.

``` r
dadaFs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 1010 sequence variants were inferred from 37907 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Dada2 a déduit 1010 séquences variantes à partir des 37907 séquences
uniques de F.

``` r
dadaRs[[1]]
```

    ## dada-class: object describing DADA2 denoising results
    ## 869 sequence variants were inferred from 45486 input unique sequences.
    ## Key parameters: OMEGA_A = 1e-40, OMEGA_C = 1e-40, BAND_SIZE = 16

Dada2 a déduit 869 séquences variantes à partir des 45486 séquences
uniques de R.

## Alignement des séquences Forward et Reverse

``` r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
```

    ## 117318 paired-reads (in 5196 unique pairings) successfully merged out of 141000 (in 21451 pairings) input.

    ## 138940 paired-reads (in 4296 unique pairings) successfully merged out of 156462 (in 15709 pairings) input.

    ## 142188 paired-reads (in 6989 unique pairings) successfully merged out of 171439 (in 27056 pairings) input.

    ## 67622 paired-reads (in 2721 unique pairings) successfully merged out of 77764 (in 9556 pairings) input.

    ## 83613 paired-reads (in 3458 unique pairings) successfully merged out of 102224 (in 16304 pairings) input.

    ## 86212 paired-reads (in 3348 unique pairings) successfully merged out of 103447 (in 14293 pairings) input.

    ## 80661 paired-reads (in 2727 unique pairings) successfully merged out of 95866 (in 12350 pairings) input.

    ## 89385 paired-reads (in 3073 unique pairings) successfully merged out of 104354 (in 12135 pairings) input.

    ## 59716 paired-reads (in 1939 unique pairings) successfully merged out of 68711 (in 7974 pairings) input.

    ## 66157 paired-reads (in 1763 unique pairings) successfully merged out of 76701 (in 8283 pairings) input.

    ## 75048 paired-reads (in 3149 unique pairings) successfully merged out of 88514 (in 12054 pairings) input.

Formation des contigues

``` r
head(mergers[[1]])
```

    ##                                                                                                                                                                                                                                                                                                                                                                                sequence
    ## 1     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 2     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 3     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTTTGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 4     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTAGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATTAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGCGAAAGCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 5     TACGAAGGGACCTAGCGTAGTTCGGAATTACTGGGCTTAAAGAGTTCGTAGGTGGTTGAAAAAGTTGGTGGTGAAATCCCAGAGCTTAACTCTGGAACTGCCATCAAAACTTTTCAGCTAGAGTATGATAGAGGAAAGCAGAATTTCTAGTGTAGAGGTGAAATTCGTAGATATTAGAAAGAATACCAATTGCGAAGGCAGCTTTCTGGATCATTACTGACACTGAGGAACGAAAGCATGGGTAGCGAAGAGGATTAGATACCCTCGTAGTCCATGCCGTAAACGATGTGTGTTAGACGTTGGAAATTTATTTTCAGTGTCGCAGGGAAACCGATAAACACACCGCCTGGGGAGTACGACCGCAAGGTT
    ## 6 TACGAGGGGTCCTAGCGTTGTCCGGATTTACTGGGCGTAAAGGGTACGTAGGCGTTTTAATAAGTTGTATGTTAAATATCTTAGCTTAACTAAGAAAGTGCATACAAAACTGTTAAGATAGAGTTTGAGAGAGGAACGCAGAATTCATGGTGGAGCGGTGACATGCGTAGATATCATGAGGAAAGTCAAATGCGAAGGCAGCCTTCTGGCTCAAAACTGACGCTGAGGTACGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTATTTGGTGCTGGGGGATTCGACCCTTTCAGTGCCGTAGCTAACGCGATAAATACTCCGCCTGGGGACTACGATCGCAAGATT
    ##   abundance forward reverse nmatch nmismatch nindel prefer accept
    ## 1      5170       1       2     29         0      0      2   TRUE
    ## 2      4129       2       1     29         0      0      2   TRUE
    ## 3      3757       3       1     29         0      0      2   TRUE
    ## 4      2481       1       1     29         0      0      2   TRUE
    ## 5      2182       2       2     29         0      0      2   TRUE
    ## 6      2132       5       9     25         0      0      1   TRUE

## Table d’observation

``` r
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
```

    ## [1]    11 19426

Construction de la table d’OTU (ASV). Sur les 11 échantillons, 19426
tables ont été créées

``` r
table(nchar(getSequences(seqtab)))
```

    ## 
    ##  352  353  362  363  364  365  366  367  368  369  370  371  372  373  374  375 
    ##    1    1    1    1    4  183   27  165  184 5608 3594 2312 2613 2738  126 1770 
    ##  376  377  378  382  386 
    ##   90    4    1    1    2

## Elimination des séquences chimériques

``` r
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
```

    ## Identified 17869 bimeras out of 19426 input sequences.

``` r
dim(seqtab.nochim)
```

    ## [1]   11 1557

Elimination des séquences chimériques

``` r
sum(seqtab.nochim)/sum(seqtab)
```

    ## [1] 0.7769154

Ces séquences chimères représentent 23% de notre jeu de données.

## Résumé des opérations effectuées précédement

``` r
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
```

    ##                             input filtered denoisedF denoisedR merged nonchim
    ## Station5_Fond1_10sept14_   159971   145448    142931    143292 117318   87962
    ## Station5_Fond1_11mars15_   175993   160423    158128    158473 138940  111552
    ## Station5_Fond2_10sept14_   197039   177018    173601    174591 142188  103668
    ## Station5_Fond2_11mars15_    87585    79989     78618     78926  67622   54711
    ## Station5_Fond3_10sept14_   117140   106150    103806    104338  83613   64259
    ## Station5_Median1_10sept14_ 116519   106745    104811    105173  86212   65559

## Téléchargement des données Silva

``` bash
wget https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
```

    ## --2020-12-19 14:08:41--  https://zenodo.org/record/3986799/files/silva_nr99_v138_train_set.fa.gz
    ## Resolving zenodo.org (zenodo.org)... 137.138.76.77
    ## Connecting to zenodo.org (zenodo.org)|137.138.76.77|:443... connected.
    ## HTTP request sent, awaiting response... 200 OK
    ## Length: 137973851 (132M) [application/octet-stream]
    ## Saving to: ‘silva_nr99_v138_train_set.fa.gz.3’
    ## 
    ##      0K .......... .......... .......... .......... ..........  0% 7.29M 18s
    ##     50K .......... .......... .......... .......... ..........  0% 3.33M 29s
    ##    100K .......... .......... .......... .......... ..........  0% 13.1M 23s
    ##    150K .......... .......... .......... .......... ..........  0% 11.9M 20s
    ##    200K .......... .......... .......... .......... ..........  0% 76.5M 16s
    ##    250K .......... .......... .......... .......... ..........  0% 95.0M 14s
    ##    300K .......... .......... .......... .......... ..........  0% 15.1M 13s
    ##    350K .......... .......... .......... .......... ..........  0% 86.4M 11s
    ##    400K .......... .......... .......... .......... ..........  0% 22.0M 11s
    ##    450K .......... .......... .......... .......... ..........  0% 52.7M 10s
    ##    500K .......... .......... .......... .......... ..........  0% 84.3M 9s
    ##    550K .......... .......... .......... .......... ..........  0% 21.5M 9s
    ##    600K .......... .......... .......... .......... ..........  0% 50.6M 8s
    ##    650K .......... .......... .......... .......... ..........  0% 78.5M 8s
    ##    700K .......... .......... .......... .......... ..........  0% 18.3M 8s
    ##    750K .......... .......... .......... .......... ..........  0% 67.1M 8s
    ##    800K .......... .......... .......... .......... ..........  0% 25.0M 7s
    ##    850K .......... .......... .......... .......... ..........  0% 74.3M 7s
    ##    900K .......... .......... .......... .......... ..........  0% 35.8M 7s
    ##    950K .......... .......... .......... .......... ..........  0% 25.8M 7s
    ##   1000K .......... .......... .......... .......... ..........  0% 29.8M 7s
    ##   1050K .......... .......... .......... .......... ..........  0% 26.0M 7s
    ##   1100K .......... .......... .......... .......... ..........  0% 69.9M 6s
    ##   1150K .......... .......... .......... .......... ..........  0% 30.5M 6s
    ##   1200K .......... .......... .......... .......... ..........  0% 26.2M 6s
    ##   1250K .......... .......... .......... .......... ..........  0% 31.1M 6s
    ##   1300K .......... .......... .......... .......... ..........  1% 30.0M 6s
    ##   1350K .......... .......... .......... .......... ..........  1% 74.7M 6s
    ##   1400K .......... .......... .......... .......... ..........  1% 26.0M 6s
    ##   1450K .......... .......... .......... .......... ..........  1% 24.6M 6s
    ##   1500K .......... .......... .......... .......... ..........  1% 28.7M 6s
    ##   1550K .......... .......... .......... .......... ..........  1% 90.1M 6s
    ##   1600K .......... .......... .......... .......... ..........  1% 23.8M 6s
    ##   1650K .......... .......... .......... .......... ..........  1% 52.2M 6s
    ##   1700K .......... .......... .......... .......... ..........  1% 18.8M 6s
    ##   1750K .......... .......... .......... .......... ..........  1% 74.0M 6s
    ##   1800K .......... .......... .......... .......... ..........  1% 76.1M 5s
    ##   1850K .......... .......... .......... .......... ..........  1% 16.0M 5s
    ##   1900K .......... .......... .......... .......... ..........  1% 71.8M 5s
    ##   1950K .......... .......... .......... .......... ..........  1% 16.6M 5s
    ##   2000K .......... .......... .......... .......... ..........  1% 68.4M 5s
    ##   2050K .......... .......... .......... .......... ..........  1% 95.3M 5s
    ##   2100K .......... .......... .......... .......... ..........  1% 14.4M 5s
    ##   2150K .......... .......... .......... .......... ..........  1% 91.1M 5s
    ##   2200K .......... .......... .......... .......... ..........  1% 17.7M 5s
    ##   2250K .......... .......... .......... .......... ..........  1% 63.2M 5s
    ##   2300K .......... .......... .......... .......... ..........  1% 69.7M 5s
    ##   2350K .......... .......... .......... .......... ..........  1% 80.2M 5s
    ##   2400K .......... .......... .......... .......... ..........  1% 21.4M 5s
    ##   2450K .......... .......... .......... .......... ..........  1% 28.3M 5s
    ##   2500K .......... .......... .......... .......... ..........  1% 71.0M 5s
    ##   2550K .......... .......... .......... .......... ..........  1% 95.7M 5s
    ##   2600K .......... .......... .......... .......... ..........  1% 17.7M 5s
    ##   2650K .......... .......... .......... .......... ..........  2% 65.9M 5s
    ##   2700K .......... .......... .......... .......... ..........  2% 83.7M 5s
    ##   2750K .......... .......... .......... .......... ..........  2%  104M 5s
    ##   2800K .......... .......... .......... .......... ..........  2% 23.2M 5s
    ##   2850K .......... .......... .......... .......... ..........  2% 85.5M 5s
    ##   2900K .......... .......... .......... .......... ..........  2% 62.7M 5s
    ##   2950K .......... .......... .......... .......... ..........  2% 36.4M 5s
    ##   3000K .......... .......... .......... .......... ..........  2% 69.1M 5s
    ##   3050K .......... .......... .......... .......... ..........  2% 22.8M 5s
    ##   3100K .......... .......... .......... .......... ..........  2% 88.7M 5s
    ##   3150K .......... .......... .......... .......... ..........  2% 41.4M 5s
    ##   3200K .......... .......... .......... .......... ..........  2% 87.7M 5s
    ##   3250K .......... .......... .......... .......... ..........  2% 28.8M 5s
    ##   3300K .......... .......... .......... .......... ..........  2% 69.5M 4s
    ##   3350K .......... .......... .......... .......... ..........  2% 68.4M 4s
    ##   3400K .......... .......... .......... .......... ..........  2% 19.4M 4s
    ##   3450K .......... .......... .......... .......... ..........  2% 78.5M 4s
    ##   3500K .......... .......... .......... .......... ..........  2% 73.9M 4s
    ##   3550K .......... .......... .......... .......... ..........  2% 55.8M 4s
    ##   3600K .......... .......... .......... .......... ..........  2% 78.9M 4s
    ##   3650K .......... .......... .......... .......... ..........  2% 37.7M 4s
    ##   3700K .......... .......... .......... .......... ..........  2% 20.4M 4s
    ##   3750K .......... .......... .......... .......... ..........  2% 99.3M 4s
    ##   3800K .......... .......... .......... .......... ..........  2% 52.9M 4s
    ##   3850K .......... .......... .......... .......... ..........  2% 20.4M 4s
    ##   3900K .......... .......... .......... .......... ..........  2% 75.7M 4s
    ##   3950K .......... .......... .......... .......... ..........  2%  104M 4s
    ##   4000K .......... .......... .......... .......... ..........  3% 67.1M 4s
    ##   4050K .......... .......... .......... .......... ..........  3% 89.8M 4s
    ##   4100K .......... .......... .......... .......... ..........  3% 3.73M 5s
    ##   4150K .......... .......... .......... .......... ..........  3%  102M 4s
    ##   4200K .......... .......... .......... .......... ..........  3%  108M 4s
    ##   4250K .......... .......... .......... .......... ..........  3%  129M 4s
    ##   4300K .......... .......... .......... .......... ..........  3% 11.3M 4s
    ##   4350K .......... .......... .......... .......... ..........  3% 71.4M 4s
    ##   4400K .......... .......... .......... .......... ..........  3% 82.7M 4s
    ##   4450K .......... .......... .......... .......... ..........  3% 50.2M 4s
    ##   4500K .......... .......... .......... .......... ..........  3% 81.0M 4s
    ##   4550K .......... .......... .......... .......... ..........  3% 29.8M 4s
    ##   4600K .......... .......... .......... .......... ..........  3% 36.0M 4s
    ##   4650K .......... .......... .......... .......... ..........  3% 78.2M 4s
    ##   4700K .......... .......... .......... .......... ..........  3% 41.9M 4s
    ##   4750K .......... .......... .......... .......... ..........  3% 38.2M 4s
    ##   4800K .......... .......... .......... .......... ..........  3% 40.3M 4s
    ##   4850K .......... .......... .......... .......... ..........  3% 54.9M 4s
    ##   4900K .......... .......... .......... .......... ..........  3% 54.6M 4s
    ##   4950K .......... .......... .......... .......... ..........  3% 42.4M 4s
    ##   5000K .......... .......... .......... .......... ..........  3% 50.8M 4s
    ##   5050K .......... .......... .......... .......... ..........  3% 56.3M 4s
    ##   5100K .......... .......... .......... .......... ..........  3% 63.8M 4s
    ##   5150K .......... .......... .......... .......... ..........  3% 46.4M 4s
    ##   5200K .......... .......... .......... .......... ..........  3% 64.7M 4s
    ##   5250K .......... .......... .......... .......... ..........  3% 49.5M 4s
    ##   5300K .......... .......... .......... .......... ..........  3% 55.1M 4s
    ##   5350K .......... .......... .......... .......... ..........  4% 59.8M 4s
    ##   5400K .......... .......... .......... .......... ..........  4% 51.3M 4s
    ##   5450K .......... .......... .......... .......... ..........  4% 27.0M 4s
    ##   5500K .......... .......... .......... .......... ..........  4% 83.5M 4s
    ##   5550K .......... .......... .......... .......... ..........  4% 72.0M 4s
    ##   5600K .......... .......... .......... .......... ..........  4% 81.6M 4s
    ##   5650K .......... .......... .......... .......... ..........  4% 9.33M 4s
    ##   5700K .......... .......... .......... .......... ..........  4% 52.0M 4s
    ##   5750K .......... .......... .......... .......... ..........  4% 37.7M 4s
    ##   5800K .......... .......... .......... .......... ..........  4%  100M 4s
    ##   5850K .......... .......... .......... .......... ..........  4% 64.6M 4s
    ##   5900K .......... .......... .......... .......... ..........  4% 94.0M 4s
    ##   5950K .......... .......... .......... .......... ..........  4% 48.4M 4s
    ##   6000K .......... .......... .......... .......... ..........  4% 46.5M 4s
    ##   6050K .......... .......... .......... .......... ..........  4% 57.7M 4s
    ##   6100K .......... .......... .......... .......... ..........  4% 90.2M 4s
    ##   6150K .......... .......... .......... .......... ..........  4% 8.13M 4s
    ##   6200K .......... .......... .......... .......... ..........  4% 31.7M 4s
    ##   6250K .......... .......... .......... .......... ..........  4% 91.3M 4s
    ##   6300K .......... .......... .......... .......... ..........  4% 54.4M 4s
    ##   6350K .......... .......... .......... .......... ..........  4% 74.1M 4s
    ##   6400K .......... .......... .......... .......... ..........  4% 78.4M 4s
    ##   6450K .......... .......... .......... .......... ..........  4% 78.8M 4s
    ##   6500K .......... .......... .......... .......... ..........  4% 38.4M 4s
    ##   6550K .......... .......... .......... .......... ..........  4% 48.3M 4s
    ##   6600K .......... .......... .......... .......... ..........  4% 73.9M 4s
    ##   6650K .......... .......... .......... .......... ..........  4% 52.0M 4s
    ##   6700K .......... .......... .......... .......... ..........  5% 82.5M 4s
    ##   6750K .......... .......... .......... .......... ..........  5% 38.3M 4s
    ##   6800K .......... .......... .......... .......... ..........  5% 34.6M 4s
    ##   6850K .......... .......... .......... .......... ..........  5% 86.8M 4s
    ##   6900K .......... .......... .......... .......... ..........  5% 75.4M 4s
    ##   6950K .......... .......... .......... .......... ..........  5% 55.1M 4s
    ##   7000K .......... .......... .......... .......... ..........  5% 44.9M 4s
    ##   7050K .......... .......... .......... .......... ..........  5% 48.0M 4s
    ##   7100K .......... .......... .......... .......... ..........  5% 77.9M 4s
    ##   7150K .......... .......... .......... .......... ..........  5% 52.9M 4s
    ##   7200K .......... .......... .......... .......... ..........  5% 47.4M 4s
    ##   7250K .......... .......... .......... .......... ..........  5% 47.7M 4s
    ##   7300K .......... .......... .......... .......... ..........  5% 30.6M 4s
    ##   7350K .......... .......... .......... .......... ..........  5% 72.0M 4s
    ##   7400K .......... .......... .......... .......... ..........  5% 41.1M 4s
    ##   7450K .......... .......... .......... .......... ..........  5% 45.8M 4s
    ##   7500K .......... .......... .......... .......... ..........  5% 34.2M 4s
    ##   7550K .......... .......... .......... .......... ..........  5% 83.4M 4s
    ##   7600K .......... .......... .......... .......... ..........  5% 37.3M 4s
    ##   7650K .......... .......... .......... .......... ..........  5% 88.6M 4s
    ##   7700K .......... .......... .......... .......... ..........  5% 66.7M 4s
    ##   7750K .......... .......... .......... .......... ..........  5% 76.6M 4s
    ##   7800K .......... .......... .......... .......... ..........  5% 90.6M 4s
    ##   7850K .......... .......... .......... .......... ..........  5% 79.5M 4s
    ##   7900K .......... .......... .......... .......... ..........  5% 82.8M 4s
    ##   7950K .......... .......... .......... .......... ..........  5% 83.7M 4s
    ##   8000K .......... .......... .......... .......... ..........  5% 26.1M 4s
    ##   8050K .......... .......... .......... .......... ..........  6% 25.8M 4s
    ##   8100K .......... .......... .......... .......... ..........  6% 54.2M 4s
    ##   8150K .......... .......... .......... .......... ..........  6%  102M 4s
    ##   8200K .......... .......... .......... .......... ..........  6% 83.2M 4s
    ##   8250K .......... .......... .......... .......... ..........  6% 90.6M 4s
    ##   8300K .......... .......... .......... .......... ..........  6% 24.5M 4s
    ##   8350K .......... .......... .......... .......... ..........  6% 94.2M 4s
    ##   8400K .......... .......... .......... .......... ..........  6% 85.8M 3s
    ##   8450K .......... .......... .......... .......... ..........  6%  100M 3s
    ##   8500K .......... .......... .......... .......... ..........  6% 83.3M 3s
    ##   8550K .......... .......... .......... .......... ..........  6% 35.3M 3s
    ##   8600K .......... .......... .......... .......... ..........  6% 53.2M 3s
    ##   8650K .......... .......... .......... .......... ..........  6% 29.4M 3s
    ##   8700K .......... .......... .......... .......... ..........  6% 67.2M 3s
    ##   8750K .......... .......... .......... .......... ..........  6% 40.6M 3s
    ##   8800K .......... .......... .......... .......... ..........  6% 43.6M 3s
    ##   8850K .......... .......... .......... .......... ..........  6% 78.4M 3s
    ##   8900K .......... .......... .......... .......... ..........  6% 42.5M 3s
    ##   8950K .......... .......... .......... .......... ..........  6% 54.9M 3s
    ##   9000K .......... .......... .......... .......... ..........  6% 50.2M 3s
    ##   9050K .......... .......... .......... .......... ..........  6% 48.4M 3s
    ##   9100K .......... .......... .......... .......... ..........  6% 42.9M 3s
    ##   9150K .......... .......... .......... .......... ..........  6% 64.3M 3s
    ##   9200K .......... .......... .......... .......... ..........  6% 42.4M 3s
    ##   9250K .......... .......... .......... .......... ..........  6% 61.5M 3s
    ##   9300K .......... .......... .......... .......... ..........  6% 80.1M 3s
    ##   9350K .......... .......... .......... .......... ..........  6% 55.3M 3s
    ##   9400K .......... .......... .......... .......... ..........  7% 70.4M 3s
    ##   9450K .......... .......... .......... .......... ..........  7% 61.8M 3s
    ##   9500K .......... .......... .......... .......... ..........  7% 77.8M 3s
    ##   9550K .......... .......... .......... .......... ..........  7% 72.5M 3s
    ##   9600K .......... .......... .......... .......... ..........  7% 75.1M 3s
    ##   9650K .......... .......... .......... .......... ..........  7% 87.0M 3s
    ##   9700K .......... .......... .......... .......... ..........  7% 71.5M 3s
    ##   9750K .......... .......... .......... .......... ..........  7% 83.0M 3s
    ##   9800K .......... .......... .......... .......... ..........  7% 20.7M 3s
    ##   9850K .......... .......... .......... .......... ..........  7% 44.1M 3s
    ##   9900K .......... .......... .......... .......... ..........  7% 82.5M 3s
    ##   9950K .......... .......... .......... .......... ..........  7% 97.9M 3s
    ##  10000K .......... .......... .......... .......... ..........  7% 81.0M 3s
    ##  10050K .......... .......... .......... .......... ..........  7% 82.9M 3s
    ##  10100K .......... .......... .......... .......... ..........  7% 19.9M 3s
    ##  10150K .......... .......... .......... .......... ..........  7% 89.8M 3s
    ##  10200K .......... .......... .......... .......... ..........  7% 72.0M 3s
    ##  10250K .......... .......... .......... .......... ..........  7% 96.1M 3s
    ##  10300K .......... .......... .......... .......... ..........  7% 45.1M 3s
    ##  10350K .......... .......... .......... .......... ..........  7% 86.7M 3s
    ##  10400K .......... .......... .......... .......... ..........  7% 50.8M 3s
    ##  10450K .......... .......... .......... .......... ..........  7%  101M 3s
    ##  10500K .......... .......... .......... .......... ..........  7% 65.9M 3s
    ##  10550K .......... .......... .......... .......... ..........  7% 67.9M 3s
    ##  10600K .......... .......... .......... .......... ..........  7% 63.6M 3s
    ##  10650K .......... .......... .......... .......... ..........  7% 79.5M 3s
    ##  10700K .......... .......... .......... .......... ..........  7% 37.8M 3s
    ##  10750K .......... .......... .......... .......... ..........  8% 30.5M 3s
    ##  10800K .......... .......... .......... .......... ..........  8% 92.1M 3s
    ##  10850K .......... .......... .......... .......... ..........  8%  104M 3s
    ##  10900K .......... .......... .......... .......... ..........  8% 72.7M 3s
    ##  10950K .......... .......... .......... .......... ..........  8% 48.4M 3s
    ##  11000K .......... .......... .......... .......... ..........  8% 57.9M 3s
    ##  11050K .......... .......... .......... .......... ..........  8% 49.6M 3s
    ##  11100K .......... .......... .......... .......... ..........  8% 49.3M 3s
    ##  11150K .......... .......... .......... .......... ..........  8% 59.4M 3s
    ##  11200K .......... .......... .......... .......... ..........  8% 89.6M 3s
    ##  11250K .......... .......... .......... .......... ..........  8% 72.4M 3s
    ##  11300K .......... .......... .......... .......... ..........  8%  106M 3s
    ##  11350K .......... .......... .......... .......... ..........  8%  107M 3s
    ##  11400K .......... .......... .......... .......... ..........  8% 42.5M 3s
    ##  11450K .......... .......... .......... .......... ..........  8% 94.9M 3s
    ##  11500K .......... .......... .......... .......... ..........  8%  108M 3s
    ##  11550K .......... .......... .......... .......... ..........  8% 36.1M 3s
    ##  11600K .......... .......... .......... .......... ..........  8% 48.9M 3s
    ##  11650K .......... .......... .......... .......... ..........  8%  117M 3s
    ##  11700K .......... .......... .......... .......... ..........  8% 24.0M 3s
    ##  11750K .......... .......... .......... .......... ..........  8% 83.9M 3s
    ##  11800K .......... .......... .......... .......... ..........  8%  120M 3s
    ##  11850K .......... .......... .......... .......... ..........  8%  107M 3s
    ##  11900K .......... .......... .......... .......... ..........  8% 72.1M 3s
    ##  11950K .......... .......... .......... .......... ..........  8% 69.7M 3s
    ##  12000K .......... .......... .......... .......... ..........  8% 79.8M 3s
    ##  12050K .......... .......... .......... .......... ..........  8% 66.9M 3s
    ##  12100K .......... .......... .......... .......... ..........  9% 28.8M 3s
    ##  12150K .......... .......... .......... .......... ..........  9% 38.0M 3s
    ##  12200K .......... .......... .......... .......... ..........  9% 98.4M 3s
    ##  12250K .......... .......... .......... .......... ..........  9%  100M 3s
    ##  12300K .......... .......... .......... .......... ..........  9%  106M 3s
    ##  12350K .......... .......... .......... .......... ..........  9% 42.5M 3s
    ##  12400K .......... .......... .......... .......... ..........  9%  119M 3s
    ##  12450K .......... .......... .......... .......... ..........  9%  133M 3s
    ##  12500K .......... .......... .......... .......... ..........  9% 63.1M 3s
    ##  12550K .......... .......... .......... .......... ..........  9% 63.8M 3s
    ##  12600K .......... .......... .......... .......... ..........  9% 66.7M 3s
    ##  12650K .......... .......... .......... .......... ..........  9% 86.1M 3s
    ##  12700K .......... .......... .......... .......... ..........  9%  105M 3s
    ##  12750K .......... .......... .......... .......... ..........  9%  147M 3s
    ##  12800K .......... .......... .......... .......... ..........  9% 50.8M 3s
    ##  12850K .......... .......... .......... .......... ..........  9% 52.7M 3s
    ##  12900K .......... .......... .......... .......... ..........  9% 72.8M 3s
    ##  12950K .......... .......... .......... .......... ..........  9% 45.2M 3s
    ##  13000K .......... .......... .......... .......... ..........  9%  115M 3s
    ##  13050K .......... .......... .......... .......... ..........  9% 23.6M 3s
    ##  13100K .......... .......... .......... .......... ..........  9% 47.1M 3s
    ##  13150K .......... .......... .......... .......... ..........  9% 60.6M 3s
    ##  13200K .......... .......... .......... .......... ..........  9% 23.1M 3s
    ##  13250K .......... .......... .......... .......... ..........  9%  138M 3s
    ##  13300K .......... .......... .......... .......... ..........  9% 56.9M 3s
    ##  13350K .......... .......... .......... .......... ..........  9% 63.2M 3s
    ##  13400K .......... .......... .......... .......... ..........  9% 99.3M 3s
    ##  13450K .......... .......... .......... .......... .......... 10%  119M 3s
    ##  13500K .......... .......... .......... .......... .......... 10% 75.5M 3s
    ##  13550K .......... .......... .......... .......... .......... 10% 77.8M 3s
    ##  13600K .......... .......... .......... .......... .......... 10% 61.6M 3s
    ##  13650K .......... .......... .......... .......... .......... 10% 74.3M 3s
    ##  13700K .......... .......... .......... .......... .......... 10% 63.2M 3s
    ##  13750K .......... .......... .......... .......... .......... 10% 87.6M 3s
    ##  13800K .......... .......... .......... .......... .......... 10% 72.8M 3s
    ##  13850K .......... .......... .......... .......... .......... 10% 78.2M 3s
    ##  13900K .......... .......... .......... .......... .......... 10% 52.0M 3s
    ##  13950K .......... .......... .......... .......... .......... 10%  146M 3s
    ##  14000K .......... .......... .......... .......... .......... 10% 87.2M 3s
    ##  14050K .......... .......... .......... .......... .......... 10% 56.2M 3s
    ##  14100K .......... .......... .......... .......... .......... 10%  112M 3s
    ##  14150K .......... .......... .......... .......... .......... 10% 17.1M 3s
    ##  14200K .......... .......... .......... .......... .......... 10% 64.9M 3s
    ##  14250K .......... .......... .......... .......... .......... 10% 18.7M 3s
    ##  14300K .......... .......... .......... .......... .......... 10% 48.1M 3s
    ##  14350K .......... .......... .......... .......... .......... 10%  127M 3s
    ##  14400K .......... .......... .......... .......... .......... 10%  135M 3s
    ##  14450K .......... .......... .......... .......... .......... 10%  153M 3s
    ##  14500K .......... .......... .......... .......... .......... 10% 83.7M 3s
    ##  14550K .......... .......... .......... .......... .......... 10%  159M 3s
    ##  14600K .......... .......... .......... .......... .......... 10% 23.2M 3s
    ##  14650K .......... .......... .......... .......... .......... 10% 78.4M 3s
    ##  14700K .......... .......... .......... .......... .......... 10% 43.8M 3s
    ##  14750K .......... .......... .......... .......... .......... 10%  147M 3s
    ##  14800K .......... .......... .......... .......... .......... 11%  103M 3s
    ##  14850K .......... .......... .......... .......... .......... 11%  159M 3s
    ##  14900K .......... .......... .......... .......... .......... 11% 31.3M 3s
    ##  14950K .......... .......... .......... .......... .......... 11% 69.2M 3s
    ##  15000K .......... .......... .......... .......... .......... 11% 74.3M 3s
    ##  15050K .......... .......... .......... .......... .......... 11%  100M 3s
    ##  15100K .......... .......... .......... .......... .......... 11%  139M 3s
    ##  15150K .......... .......... .......... .......... .......... 11%  138M 3s
    ##  15200K .......... .......... .......... .......... .......... 11% 15.8M 3s
    ##  15250K .......... .......... .......... .......... .......... 11%  123M 3s
    ##  15300K .......... .......... .......... .......... .......... 11% 83.1M 3s
    ##  15350K .......... .......... .......... .......... .......... 11%  166M 3s
    ##  15400K .......... .......... .......... .......... .......... 11%  144M 3s
    ##  15450K .......... .......... .......... .......... .......... 11%  168M 3s
    ##  15500K .......... .......... .......... .......... .......... 11% 23.3M 3s
    ##  15550K .......... .......... .......... .......... .......... 11%  109M 3s
    ##  15600K .......... .......... .......... .......... .......... 11% 77.3M 3s
    ##  15650K .......... .......... .......... .......... .......... 11% 64.5M 3s
    ##  15700K .......... .......... .......... .......... .......... 11%  145M 3s
    ##  15750K .......... .......... .......... .......... .......... 11%  180M 3s
    ##  15800K .......... .......... .......... .......... .......... 11% 10.4M 3s
    ##  15850K .......... .......... .......... .......... .......... 11% 54.1M 3s
    ##  15900K .......... .......... .......... .......... .......... 11% 59.7M 3s
    ##  15950K .......... .......... .......... .......... .......... 11% 42.4M 3s
    ##  16000K .......... .......... .......... .......... .......... 11%  148M 3s
    ##  16050K .......... .......... .......... .......... .......... 11%  181M 3s
    ##  16100K .......... .......... .......... .......... .......... 11% 7.38M 3s
    ##  16150K .......... .......... .......... .......... .......... 12% 32.3M 3s
    ##  16200K .......... .......... .......... .......... .......... 12% 50.4M 3s
    ##  16250K .......... .......... .......... .......... .......... 12% 52.6M 3s
    ##  16300K .......... .......... .......... .......... .......... 12% 98.1M 3s
    ##  16350K .......... .......... .......... .......... .......... 12% 53.3M 3s
    ##  16400K .......... .......... .......... .......... .......... 12%  130M 3s
    ##  16450K .......... .......... .......... .......... .......... 12%  145M 3s
    ##  16500K .......... .......... .......... .......... .......... 12% 74.3M 3s
    ##  16550K .......... .......... .......... .......... .......... 12% 53.9M 3s
    ##  16600K .......... .......... .......... .......... .......... 12% 52.4M 3s
    ##  16650K .......... .......... .......... .......... .......... 12% 92.3M 3s
    ##  16700K .......... .......... .......... .......... .......... 12% 43.3M 3s
    ##  16750K .......... .......... .......... .......... .......... 12% 87.3M 3s
    ##  16800K .......... .......... .......... .......... .......... 12% 45.6M 3s
    ##  16850K .......... .......... .......... .......... .......... 12% 43.3M 3s
    ##  16900K .......... .......... .......... .......... .......... 12%  147M 3s
    ##  16950K .......... .......... .......... .......... .......... 12%  115M 3s
    ##  17000K .......... .......... .......... .......... .......... 12%  150M 3s
    ##  17050K .......... .......... .......... .......... .......... 12% 81.6M 3s
    ##  17100K .......... .......... .......... .......... .......... 12% 21.9M 3s
    ##  17150K .......... .......... .......... .......... .......... 12%  155M 3s
    ##  17200K .......... .......... .......... .......... .......... 12%  158M 3s
    ##  17250K .......... .......... .......... .......... .......... 12%  142M 3s
    ##  17300K .......... .......... .......... .......... .......... 12% 85.8M 3s
    ##  17350K .......... .......... .......... .......... .......... 12% 93.9M 3s
    ##  17400K .......... .......... .......... .......... .......... 12%  141M 3s
    ##  17450K .......... .......... .......... .......... .......... 12% 33.6M 3s
    ##  17500K .......... .......... .......... .......... .......... 13% 69.2M 3s
    ##  17550K .......... .......... .......... .......... .......... 13% 40.6M 3s
    ##  17600K .......... .......... .......... .......... .......... 13%  148M 3s
    ##  17650K .......... .......... .......... .......... .......... 13% 38.3M 3s
    ##  17700K .......... .......... .......... .......... .......... 13% 69.1M 3s
    ##  17750K .......... .......... .......... .......... .......... 13%  143M 3s
    ##  17800K .......... .......... .......... .......... .......... 13% 91.1M 3s
    ##  17850K .......... .......... .......... .......... .......... 13%  110M 3s
    ##  17900K .......... .......... .......... .......... .......... 13% 80.9M 3s
    ##  17950K .......... .......... .......... .......... .......... 13% 55.6M 3s
    ##  18000K .......... .......... .......... .......... .......... 13%  105M 3s
    ##  18050K .......... .......... .......... .......... .......... 13% 85.3M 3s
    ##  18100K .......... .......... .......... .......... .......... 13% 45.6M 3s
    ##  18150K .......... .......... .......... .......... .......... 13%  114M 3s
    ##  18200K .......... .......... .......... .......... .......... 13% 65.1M 3s
    ##  18250K .......... .......... .......... .......... .......... 13% 75.8M 3s
    ##  18300K .......... .......... .......... .......... .......... 13%  134M 3s
    ##  18350K .......... .......... .......... .......... .......... 13% 63.4M 3s
    ##  18400K .......... .......... .......... .......... .......... 13% 69.3M 3s
    ##  18450K .......... .......... .......... .......... .......... 13% 41.7M 3s
    ##  18500K .......... .......... .......... .......... .......... 13% 69.2M 3s
    ##  18550K .......... .......... .......... .......... .......... 13% 51.1M 3s
    ##  18600K .......... .......... .......... .......... .......... 13%  113M 3s
    ##  18650K .......... .......... .......... .......... .......... 13%  152M 3s
    ##  18700K .......... .......... .......... .......... .......... 13% 87.1M 3s
    ##  18750K .......... .......... .......... .......... .......... 13% 58.7M 3s
    ##  18800K .......... .......... .......... .......... .......... 13% 89.5M 3s
    ##  18850K .......... .......... .......... .......... .......... 14% 7.93M 3s
    ##  18900K .......... .......... .......... .......... .......... 14% 30.3M 3s
    ##  18950K .......... .......... .......... .......... .......... 14%  139M 3s
    ##  19000K .......... .......... .......... .......... .......... 14% 64.6M 3s
    ##  19050K .......... .......... .......... .......... .......... 14% 59.2M 3s
    ##  19100K .......... .......... .......... .......... .......... 14%  124M 3s
    ##  19150K .......... .......... .......... .......... .......... 14%  145M 3s
    ##  19200K .......... .......... .......... .......... .......... 14%  105M 3s
    ##  19250K .......... .......... .......... .......... .......... 14% 86.2M 3s
    ##  19300K .......... .......... .......... .......... .......... 14% 78.2M 3s
    ##  19350K .......... .......... .......... .......... .......... 14% 35.3M 3s
    ##  19400K .......... .......... .......... .......... .......... 14% 91.0M 3s
    ##  19450K .......... .......... .......... .......... .......... 14%  166M 2s
    ##  19500K .......... .......... .......... .......... .......... 14% 7.50M 3s
    ##  19550K .......... .......... .......... .......... .......... 14%  117M 3s
    ##  19600K .......... .......... .......... .......... .......... 14% 73.2M 3s
    ##  19650K .......... .......... .......... .......... .......... 14% 48.6M 3s
    ##  19700K .......... .......... .......... .......... .......... 14%  100M 3s
    ##  19750K .......... .......... .......... .......... .......... 14% 73.4M 3s
    ##  19800K .......... .......... .......... .......... .......... 14% 76.8M 3s
    ##  19850K .......... .......... .......... .......... .......... 14%  133M 2s
    ##  19900K .......... .......... .......... .......... .......... 14% 65.1M 2s
    ##  19950K .......... .......... .......... .......... .......... 14%  141M 2s
    ##  20000K .......... .......... .......... .......... .......... 14% 15.5M 3s
    ##  20050K .......... .......... .......... .......... .......... 14%  162M 2s
    ##  20100K .......... .......... .......... .......... .......... 14%  123M 2s
    ##  20150K .......... .......... .......... .......... .......... 14%  116M 2s
    ##  20200K .......... .......... .......... .......... .......... 15% 88.0M 2s
    ##  20250K .......... .......... .......... .......... .......... 15% 20.5M 2s
    ##  20300K .......... .......... .......... .......... .......... 15% 42.9M 2s
    ##  20350K .......... .......... .......... .......... .......... 15% 50.2M 2s
    ##  20400K .......... .......... .......... .......... .......... 15%  106M 2s
    ##  20450K .......... .......... .......... .......... .......... 15%  156M 2s
    ##  20500K .......... .......... .......... .......... .......... 15%  150M 2s
    ##  20550K .......... .......... .......... .......... .......... 15%  154M 2s
    ##  20600K .......... .......... .......... .......... .......... 15% 67.8M 2s
    ##  20650K .......... .......... .......... .......... .......... 15% 80.0M 2s
    ##  20700K .......... .......... .......... .......... .......... 15%  147M 2s
    ##  20750K .......... .......... .......... .......... .......... 15%  175M 2s
    ##  20800K .......... .......... .......... .......... .......... 15%  140M 2s
    ##  20850K .......... .......... .......... .......... .......... 15% 65.7M 2s
    ##  20900K .......... .......... .......... .......... .......... 15%  136M 2s
    ##  20950K .......... .......... .......... .......... .......... 15%  121M 2s
    ##  21000K .......... .......... .......... .......... .......... 15% 79.8M 2s
    ##  21050K .......... .......... .......... .......... .......... 15% 58.0M 2s
    ##  21100K .......... .......... .......... .......... .......... 15% 81.0M 2s
    ##  21150K .......... .......... .......... .......... .......... 15% 77.0M 2s
    ##  21200K .......... .......... .......... .......... .......... 15% 42.9M 2s
    ##  21250K .......... .......... .......... .......... .......... 15%  155M 2s
    ##  21300K .......... .......... .......... .......... .......... 15% 66.7M 2s
    ##  21350K .......... .......... .......... .......... .......... 15% 56.9M 2s
    ##  21400K .......... .......... .......... .......... .......... 15% 96.3M 2s
    ##  21450K .......... .......... .......... .......... .......... 15% 83.2M 2s
    ##  21500K .......... .......... .......... .......... .......... 15% 61.7M 2s
    ##  21550K .......... .......... .......... .......... .......... 16% 98.1M 2s
    ##  21600K .......... .......... .......... .......... .......... 16% 96.7M 2s
    ##  21650K .......... .......... .......... .......... .......... 16%  114M 2s
    ##  21700K .......... .......... .......... .......... .......... 16% 62.1M 2s
    ##  21750K .......... .......... .......... .......... .......... 16% 96.0M 2s
    ##  21800K .......... .......... .......... .......... .......... 16% 63.6M 2s
    ##  21850K .......... .......... .......... .......... .......... 16% 60.0M 2s
    ##  21900K .......... .......... .......... .......... .......... 16%  120M 2s
    ##  21950K .......... .......... .......... .......... .......... 16%  139M 2s
    ##  22000K .......... .......... .......... .......... .......... 16% 50.1M 2s
    ##  22050K .......... .......... .......... .......... .......... 16% 98.2M 2s
    ##  22100K .......... .......... .......... .......... .......... 16% 42.7M 2s
    ##  22150K .......... .......... .......... .......... .......... 16% 74.7M 2s
    ##  22200K .......... .......... .......... .......... .......... 16% 81.2M 2s
    ##  22250K .......... .......... .......... .......... .......... 16% 90.0M 2s
    ##  22300K .......... .......... .......... .......... .......... 16% 43.6M 2s
    ##  22350K .......... .......... .......... .......... .......... 16% 45.9M 2s
    ##  22400K .......... .......... .......... .......... .......... 16% 46.9M 2s
    ##  22450K .......... .......... .......... .......... .......... 16% 39.7M 2s
    ##  22500K .......... .......... .......... .......... .......... 16% 46.9M 2s
    ##  22550K .......... .......... .......... .......... .......... 16% 69.9M 2s
    ##  22600K .......... .......... .......... .......... .......... 16% 75.7M 2s
    ##  22650K .......... .......... .......... .......... .......... 16% 54.7M 2s
    ##  22700K .......... .......... .......... .......... .......... 16% 39.9M 2s
    ##  22750K .......... .......... .......... .......... .......... 16% 32.6M 2s
    ##  22800K .......... .......... .......... .......... .......... 16% 57.2M 2s
    ##  22850K .......... .......... .......... .......... .......... 16% 63.4M 2s
    ##  22900K .......... .......... .......... .......... .......... 17% 58.2M 2s
    ##  22950K .......... .......... .......... .......... .......... 17% 45.5M 2s
    ##  23000K .......... .......... .......... .......... .......... 17% 58.8M 2s
    ##  23050K .......... .......... .......... .......... .......... 17% 58.0M 2s
    ##  23100K .......... .......... .......... .......... .......... 17% 57.6M 2s
    ##  23150K .......... .......... .......... .......... .......... 17% 92.9M 2s
    ##  23200K .......... .......... .......... .......... .......... 17% 42.2M 2s
    ##  23250K .......... .......... .......... .......... .......... 17% 92.3M 2s
    ##  23300K .......... .......... .......... .......... .......... 17% 76.5M 2s
    ##  23350K .......... .......... .......... .......... .......... 17% 49.9M 2s
    ##  23400K .......... .......... .......... .......... .......... 17% 84.7M 2s
    ##  23450K .......... .......... .......... .......... .......... 17% 56.0M 2s
    ##  23500K .......... .......... .......... .......... .......... 17% 25.6M 2s
    ##  23550K .......... .......... .......... .......... .......... 17%  110M 2s
    ##  23600K .......... .......... .......... .......... .......... 17% 73.9M 2s
    ##  23650K .......... .......... .......... .......... .......... 17% 92.2M 2s
    ##  23700K .......... .......... .......... .......... .......... 17% 45.3M 2s
    ##  23750K .......... .......... .......... .......... .......... 17% 50.6M 2s
    ##  23800K .......... .......... .......... .......... .......... 17% 99.9M 2s
    ##  23850K .......... .......... .......... .......... .......... 17% 97.5M 2s
    ##  23900K .......... .......... .......... .......... .......... 17% 74.8M 2s
    ##  23950K .......... .......... .......... .......... .......... 17% 59.5M 2s
    ##  24000K .......... .......... .......... .......... .......... 17% 38.2M 2s
    ##  24050K .......... .......... .......... .......... .......... 17% 80.9M 2s
    ##  24100K .......... .......... .......... .......... .......... 17%  105M 2s
    ##  24150K .......... .......... .......... .......... .......... 17% 45.7M 2s
    ##  24200K .......... .......... .......... .......... .......... 17% 49.3M 2s
    ##  24250K .......... .......... .......... .......... .......... 18%  111M 2s
    ##  24300K .......... .......... .......... .......... .......... 18%  102M 2s
    ##  24350K .......... .......... .......... .......... .......... 18% 59.2M 2s
    ##  24400K .......... .......... .......... .......... .......... 18% 70.8M 2s
    ##  24450K .......... .......... .......... .......... .......... 18% 75.9M 2s
    ##  24500K .......... .......... .......... .......... .......... 18% 37.9M 2s
    ##  24550K .......... .......... .......... .......... .......... 18% 87.9M 2s
    ##  24600K .......... .......... .......... .......... .......... 18% 64.3M 2s
    ##  24650K .......... .......... .......... .......... .......... 18% 99.6M 2s
    ##  24700K .......... .......... .......... .......... .......... 18%  106M 2s
    ##  24750K .......... .......... .......... .......... .......... 18%  107M 2s
    ##  24800K .......... .......... .......... .......... .......... 18% 99.5M 2s
    ##  24850K .......... .......... .......... .......... .......... 18% 94.3M 2s
    ##  24900K .......... .......... .......... .......... .......... 18% 54.4M 2s
    ##  24950K .......... .......... .......... .......... .......... 18% 81.3M 2s
    ##  25000K .......... .......... .......... .......... .......... 18%  102M 2s
    ##  25050K .......... .......... .......... .......... .......... 18%  119M 2s
    ##  25100K .......... .......... .......... .......... .......... 18%  104M 2s
    ##  25150K .......... .......... .......... .......... .......... 18% 98.0M 2s
    ##  25200K .......... .......... .......... .......... .......... 18%  109M 2s
    ##  25250K .......... .......... .......... .......... .......... 18%  112M 2s
    ##  25300K .......... .......... .......... .......... .......... 18% 94.9M 2s
    ##  25350K .......... .......... .......... .......... .......... 18% 82.9M 2s
    ##  25400K .......... .......... .......... .......... .......... 18% 71.6M 2s
    ##  25450K .......... .......... .......... .......... .......... 18% 71.9M 2s
    ##  25500K .......... .......... .......... .......... .......... 18% 74.1M 2s
    ##  25550K .......... .......... .......... .......... .......... 18% 28.1M 2s
    ##  25600K .......... .......... .......... .......... .......... 19% 55.8M 2s
    ##  25650K .......... .......... .......... .......... .......... 19% 99.7M 2s
    ##  25700K .......... .......... .......... .......... .......... 19%  112M 2s
    ##  25750K .......... .......... .......... .......... .......... 19%  101M 2s
    ##  25800K .......... .......... .......... .......... .......... 19% 84.5M 2s
    ##  25850K .......... .......... .......... .......... .......... 19% 75.9M 2s
    ##  25900K .......... .......... .......... .......... .......... 19% 66.4M 2s
    ##  25950K .......... .......... .......... .......... .......... 19% 89.0M 2s
    ##  26000K .......... .......... .......... .......... .......... 19% 38.6M 2s
    ##  26050K .......... .......... .......... .......... .......... 19% 52.2M 2s
    ##  26100K .......... .......... .......... .......... .......... 19% 11.4M 2s
    ##  26150K .......... .......... .......... .......... .......... 19%  136M 2s
    ##  26200K .......... .......... .......... .......... .......... 19%  124M 2s
    ##  26250K .......... .......... .......... .......... .......... 19%  131M 2s
    ##  26300K .......... .......... .......... .......... .......... 19% 15.8M 2s
    ##  26350K .......... .......... .......... .......... .......... 19% 31.9M 2s
    ##  26400K .......... .......... .......... .......... .......... 19% 72.5M 2s
    ##  26450K .......... .......... .......... .......... .......... 19% 62.2M 2s
    ##  26500K .......... .......... .......... .......... .......... 19% 40.5M 2s
    ##  26550K .......... .......... .......... .......... .......... 19% 54.0M 2s
    ##  26600K .......... .......... .......... .......... .......... 19% 37.1M 2s
    ##  26650K .......... .......... .......... .......... .......... 19% 93.7M 2s
    ##  26700K .......... .......... .......... .......... .......... 19% 55.9M 2s
    ##  26750K .......... .......... .......... .......... .......... 19% 55.3M 2s
    ##  26800K .......... .......... .......... .......... .......... 19% 51.2M 2s
    ##  26850K .......... .......... .......... .......... .......... 19% 93.9M 2s
    ##  26900K .......... .......... .......... .......... .......... 20% 74.9M 2s
    ##  26950K .......... .......... .......... .......... .......... 20% 48.4M 2s
    ##  27000K .......... .......... .......... .......... .......... 20% 58.8M 2s
    ##  27050K .......... .......... .......... .......... .......... 20% 40.1M 2s
    ##  27100K .......... .......... .......... .......... .......... 20% 97.3M 2s
    ##  27150K .......... .......... .......... .......... .......... 20% 60.3M 2s
    ##  27200K .......... .......... .......... .......... .......... 20% 83.8M 2s
    ##  27250K .......... .......... .......... .......... .......... 20% 70.1M 2s
    ##  27300K .......... .......... .......... .......... .......... 20% 66.0M 2s
    ##  27350K .......... .......... .......... .......... .......... 20% 93.8M 2s
    ##  27400K .......... .......... .......... .......... .......... 20% 65.2M 2s
    ##  27450K .......... .......... .......... .......... .......... 20% 98.0M 2s
    ##  27500K .......... .......... .......... .......... .......... 20% 77.7M 2s
    ##  27550K .......... .......... .......... .......... .......... 20% 51.4M 2s
    ##  27600K .......... .......... .......... .......... .......... 20% 33.4M 2s
    ##  27650K .......... .......... .......... .......... .......... 20% 98.3M 2s
    ##  27700K .......... .......... .......... .......... .......... 20% 98.6M 2s
    ##  27750K .......... .......... .......... .......... .......... 20% 67.2M 2s
    ##  27800K .......... .......... .......... .......... .......... 20% 99.3M 2s
    ##  27850K .......... .......... .......... .......... .......... 20% 59.1M 2s
    ##  27900K .......... .......... .......... .......... .......... 20% 58.9M 2s
    ##  27950K .......... .......... .......... .......... .......... 20% 86.8M 2s
    ##  28000K .......... .......... .......... .......... .......... 20% 53.7M 2s
    ##  28050K .......... .......... .......... .......... .......... 20% 66.2M 2s
    ##  28100K .......... .......... .......... .......... .......... 20% 95.6M 2s
    ##  28150K .......... .......... .......... .......... .......... 20% 97.9M 2s
    ##  28200K .......... .......... .......... .......... .......... 20% 92.2M 2s
    ##  28250K .......... .......... .......... .......... .......... 21% 71.9M 2s
    ##  28300K .......... .......... .......... .......... .......... 21% 79.6M 2s
    ##  28350K .......... .......... .......... .......... .......... 21% 82.6M 2s
    ##  28400K .......... .......... .......... .......... .......... 21%  101M 2s
    ##  28450K .......... .......... .......... .......... .......... 21% 90.7M 2s
    ##  28500K .......... .......... .......... .......... .......... 21% 89.9M 2s
    ##  28550K .......... .......... .......... .......... .......... 21% 60.7M 2s
    ##  28600K .......... .......... .......... .......... .......... 21% 47.8M 2s
    ##  28650K .......... .......... .......... .......... .......... 21% 88.0M 2s
    ##  28700K .......... .......... .......... .......... .......... 21% 69.4M 2s
    ##  28750K .......... .......... .......... .......... .......... 21% 54.8M 2s
    ##  28800K .......... .......... .......... .......... .......... 21% 71.4M 2s
    ##  28850K .......... .......... .......... .......... .......... 21%  110M 2s
    ##  28900K .......... .......... .......... .......... .......... 21%  102M 2s
    ##  28950K .......... .......... .......... .......... .......... 21%  116M 2s
    ##  29000K .......... .......... .......... .......... .......... 21% 96.8M 2s
    ##  29050K .......... .......... .......... .......... .......... 21% 56.5M 2s
    ##  29100K .......... .......... .......... .......... .......... 21% 64.1M 2s
    ##  29150K .......... .......... .......... .......... .......... 21%  114M 2s
    ##  29200K .......... .......... .......... .......... .......... 21% 94.9M 2s
    ##  29250K .......... .......... .......... .......... .......... 21%  108M 2s
    ##  29300K .......... .......... .......... .......... .......... 21% 68.6M 2s
    ##  29350K .......... .......... .......... .......... .......... 21%  104M 2s
    ##  29400K .......... .......... .......... .......... .......... 21% 90.7M 2s
    ##  29450K .......... .......... .......... .......... .......... 21% 74.6M 2s
    ##  29500K .......... .......... .......... .......... .......... 21%  100M 2s
    ##  29550K .......... .......... .......... .......... .......... 21% 45.5M 2s
    ##  29600K .......... .......... .......... .......... .......... 22% 55.9M 2s
    ##  29650K .......... .......... .......... .......... .......... 22% 77.8M 2s
    ##  29700K .......... .......... .......... .......... .......... 22% 96.5M 2s
    ##  29750K .......... .......... .......... .......... .......... 22% 58.4M 2s
    ##  29800K .......... .......... .......... .......... .......... 22% 77.4M 2s
    ##  29850K .......... .......... .......... .......... .......... 22% 74.1M 2s
    ##  29900K .......... .......... .......... .......... .......... 22% 57.0M 2s
    ##  29950K .......... .......... .......... .......... .......... 22%  121M 2s
    ##  30000K .......... .......... .......... .......... .......... 22%  108M 2s
    ##  30050K .......... .......... .......... .......... .......... 22% 93.6M 2s
    ##  30100K .......... .......... .......... .......... .......... 22% 77.4M 2s
    ##  30150K .......... .......... .......... .......... .......... 22% 76.0M 2s
    ##  30200K .......... .......... .......... .......... .......... 22% 69.1M 2s
    ##  30250K .......... .......... .......... .......... .......... 22% 71.5M 2s
    ##  30300K .......... .......... .......... .......... .......... 22%  122M 2s
    ##  30350K .......... .......... .......... .......... .......... 22% 48.7M 2s
    ##  30400K .......... .......... .......... .......... .......... 22% 99.6M 2s
    ##  30450K .......... .......... .......... .......... .......... 22% 87.8M 2s
    ##  30500K .......... .......... .......... .......... .......... 22%  113M 2s
    ##  30550K .......... .......... .......... .......... .......... 22% 95.4M 2s
    ##  30600K .......... .......... .......... .......... .......... 22% 43.7M 2s
    ##  30650K .......... .......... .......... .......... .......... 22% 92.9M 2s
    ##  30700K .......... .......... .......... .......... .......... 22% 78.9M 2s
    ##  30750K .......... .......... .......... .......... .......... 22% 63.8M 2s
    ##  30800K .......... .......... .......... .......... .......... 22% 61.0M 2s
    ##  30850K .......... .......... .......... .......... .......... 22% 59.9M 2s
    ##  30900K .......... .......... .......... .......... .......... 22% 63.6M 2s
    ##  30950K .......... .......... .......... .......... .......... 23% 77.6M 2s
    ##  31000K .......... .......... .......... .......... .......... 23%  111M 2s
    ##  31050K .......... .......... .......... .......... .......... 23%  143M 2s
    ##  31100K .......... .......... .......... .......... .......... 23%  113M 2s
    ##  31150K .......... .......... .......... .......... .......... 23% 66.2M 2s
    ##  31200K .......... .......... .......... .......... .......... 23% 75.6M 2s
    ##  31250K .......... .......... .......... .......... .......... 23% 74.1M 2s
    ##  31300K .......... .......... .......... .......... .......... 23% 84.6M 2s
    ##  31350K .......... .......... .......... .......... .......... 23%  149M 2s
    ##  31400K .......... .......... .......... .......... .......... 23%  104M 2s
    ##  31450K .......... .......... .......... .......... .......... 23%  105M 2s
    ##  31500K .......... .......... .......... .......... .......... 23% 70.5M 2s
    ##  31550K .......... .......... .......... .......... .......... 23% 68.1M 2s
    ##  31600K .......... .......... .......... .......... .......... 23% 84.4M 2s
    ##  31650K .......... .......... .......... .......... .......... 23% 85.1M 2s
    ##  31700K .......... .......... .......... .......... .......... 23% 61.6M 2s
    ##  31750K .......... .......... .......... .......... .......... 23% 39.1M 2s
    ##  31800K .......... .......... .......... .......... .......... 23%  108M 2s
    ##  31850K .......... .......... .......... .......... .......... 23%  107M 2s
    ##  31900K .......... .......... .......... .......... .......... 23% 98.9M 2s
    ##  31950K .......... .......... .......... .......... .......... 23% 87.9M 2s
    ##  32000K .......... .......... .......... .......... .......... 23% 94.5M 2s
    ##  32050K .......... .......... .......... .......... .......... 23% 95.7M 2s
    ##  32100K .......... .......... .......... .......... .......... 23%  109M 2s
    ##  32150K .......... .......... .......... .......... .......... 23%  101M 2s
    ##  32200K .......... .......... .......... .......... .......... 23%  140M 2s
    ##  32250K .......... .......... .......... .......... .......... 23% 78.6M 2s
    ##  32300K .......... .......... .......... .......... .......... 24%  122M 2s
    ##  32350K .......... .......... .......... .......... .......... 24%  141M 2s
    ##  32400K .......... .......... .......... .......... .......... 24%  106M 2s
    ##  32450K .......... .......... .......... .......... .......... 24% 85.8M 2s
    ##  32500K .......... .......... .......... .......... .......... 24%  140M 2s
    ##  32550K .......... .......... .......... .......... .......... 24% 47.0M 2s
    ##  32600K .......... .......... .......... .......... .......... 24% 66.6M 2s
    ##  32650K .......... .......... .......... .......... .......... 24% 74.4M 2s
    ##  32700K .......... .......... .......... .......... .......... 24% 82.4M 2s
    ##  32750K .......... .......... .......... .......... .......... 24% 48.8M 2s
    ##  32800K .......... .......... .......... .......... .......... 24% 92.9M 2s
    ##  32850K .......... .......... .......... .......... .......... 24% 80.4M 2s
    ##  32900K .......... .......... .......... .......... .......... 24% 63.2M 2s
    ##  32950K .......... .......... .......... .......... .......... 24% 73.6M 2s
    ##  33000K .......... .......... .......... .......... .......... 24% 68.2M 2s
    ##  33050K .......... .......... .......... .......... .......... 24%  136M 2s
    ##  33100K .......... .......... .......... .......... .......... 24% 60.8M 2s
    ##  33150K .......... .......... .......... .......... .......... 24% 69.3M 2s
    ##  33200K .......... .......... .......... .......... .......... 24% 90.5M 2s
    ##  33250K .......... .......... .......... .......... .......... 24%  104M 2s
    ##  33300K .......... .......... .......... .......... .......... 24% 52.7M 2s
    ##  33350K .......... .......... .......... .......... .......... 24% 60.6M 2s
    ##  33400K .......... .......... .......... .......... .......... 24% 87.9M 2s
    ##  33450K .......... .......... .......... .......... .......... 24% 95.3M 2s
    ##  33500K .......... .......... .......... .......... .......... 24% 61.6M 2s
    ##  33550K .......... .......... .......... .......... .......... 24% 58.2M 2s
    ##  33600K .......... .......... .......... .......... .......... 24% 63.8M 2s
    ##  33650K .......... .......... .......... .......... .......... 25%  108M 2s
    ##  33700K .......... .......... .......... .......... .......... 25% 94.4M 2s
    ##  33750K .......... .......... .......... .......... .......... 25% 67.6M 2s
    ##  33800K .......... .......... .......... .......... .......... 25% 69.4M 2s
    ##  33850K .......... .......... .......... .......... .......... 25%  105M 2s
    ##  33900K .......... .......... .......... .......... .......... 25% 80.5M 2s
    ##  33950K .......... .......... .......... .......... .......... 25% 68.6M 2s
    ##  34000K .......... .......... .......... .......... .......... 25% 82.1M 2s
    ##  34050K .......... .......... .......... .......... .......... 25%  110M 2s
    ##  34100K .......... .......... .......... .......... .......... 25% 65.7M 2s
    ##  34150K .......... .......... .......... .......... .......... 25% 65.4M 2s
    ##  34200K .......... .......... .......... .......... .......... 25% 60.9M 2s
    ##  34250K .......... .......... .......... .......... .......... 25% 86.1M 2s
    ##  34300K .......... .......... .......... .......... .......... 25% 90.4M 2s
    ##  34350K .......... .......... .......... .......... .......... 25% 95.1M 2s
    ##  34400K .......... .......... .......... .......... .......... 25% 82.0M 2s
    ##  34450K .......... .......... .......... .......... .......... 25% 54.2M 2s
    ##  34500K .......... .......... .......... .......... .......... 25% 76.4M 2s
    ##  34550K .......... .......... .......... .......... .......... 25%  106M 2s
    ##  34600K .......... .......... .......... .......... .......... 25% 72.3M 2s
    ##  34650K .......... .......... .......... .......... .......... 25% 73.8M 2s
    ##  34700K .......... .......... .......... .......... .......... 25% 70.1M 2s
    ##  34750K .......... .......... .......... .......... .......... 25% 78.7M 2s
    ##  34800K .......... .......... .......... .......... .......... 25%  128M 2s
    ##  34850K .......... .......... .......... .......... .......... 25% 73.7M 2s
    ##  34900K .......... .......... .......... .......... .......... 25% 75.5M 2s
    ##  34950K .......... .......... .......... .......... .......... 25% 71.4M 2s
    ##  35000K .......... .......... .......... .......... .......... 26%  146M 2s
    ##  35050K .......... .......... .......... .......... .......... 26% 78.9M 2s
    ##  35100K .......... .......... .......... .......... .......... 26% 75.2M 2s
    ##  35150K .......... .......... .......... .......... .......... 26% 73.4M 2s
    ##  35200K .......... .......... .......... .......... .......... 26%  157M 2s
    ##  35250K .......... .......... .......... .......... .......... 26% 70.5M 2s
    ##  35300K .......... .......... .......... .......... .......... 26% 69.8M 2s
    ##  35350K .......... .......... .......... .......... .......... 26% 92.2M 2s
    ##  35400K .......... .......... .......... .......... .......... 26% 65.5M 2s
    ##  35450K .......... .......... .......... .......... .......... 26%  112M 2s
    ##  35500K .......... .......... .......... .......... .......... 26% 89.5M 2s
    ##  35550K .......... .......... .......... .......... .......... 26% 75.1M 2s
    ##  35600K .......... .......... .......... .......... .......... 26% 62.2M 2s
    ##  35650K .......... .......... .......... .......... .......... 26% 85.1M 2s
    ##  35700K .......... .......... .......... .......... .......... 26% 70.2M 2s
    ##  35750K .......... .......... .......... .......... .......... 26% 75.7M 2s
    ##  35800K .......... .......... .......... .......... .......... 26%  153M 2s
    ##  35850K .......... .......... .......... .......... .......... 26% 74.3M 2s
    ##  35900K .......... .......... .......... .......... .......... 26% 70.7M 2s
    ##  35950K .......... .......... .......... .......... .......... 26% 81.0M 2s
    ##  36000K .......... .......... .......... .......... .......... 26% 68.2M 2s
    ##  36050K .......... .......... .......... .......... .......... 26%  149M 2s
    ##  36100K .......... .......... .......... .......... .......... 26%  133M 2s
    ##  36150K .......... .......... .......... .......... .......... 26%  140M 2s
    ##  36200K .......... .......... .......... .......... .......... 26%  107M 2s
    ##  36250K .......... .......... .......... .......... .......... 26% 98.3M 2s
    ##  36300K .......... .......... .......... .......... .......... 26%  148M 2s
    ##  36350K .......... .......... .......... .......... .......... 27%  165M 2s
    ##  36400K .......... .......... .......... .......... .......... 27%  139M 2s
    ##  36450K .......... .......... .......... .......... .......... 27%  113M 2s
    ##  36500K .......... .......... .......... .......... .......... 27%  139M 2s
    ##  36550K .......... .......... .......... .......... .......... 27%  139M 2s
    ##  36600K .......... .......... .......... .......... .......... 27% 91.1M 2s
    ##  36650K .......... .......... .......... .......... .......... 27% 98.5M 2s
    ##  36700K .......... .......... .......... .......... .......... 27% 44.3M 2s
    ##  36750K .......... .......... .......... .......... .......... 27%  162M 2s
    ##  36800K .......... .......... .......... .......... .......... 27%  116M 2s
    ##  36850K .......... .......... .......... .......... .......... 27% 64.5M 2s
    ##  36900K .......... .......... .......... .......... .......... 27%  134M 2s
    ##  36950K .......... .......... .......... .......... .......... 27%  151M 2s
    ##  37000K .......... .......... .......... .......... .......... 27%  151M 2s
    ##  37050K .......... .......... .......... .......... .......... 27%  154M 2s
    ##  37100K .......... .......... .......... .......... .......... 27% 55.4M 2s
    ##  37150K .......... .......... .......... .......... .......... 27%  161M 2s
    ##  37200K .......... .......... .......... .......... .......... 27% 76.0M 2s
    ##  37250K .......... .......... .......... .......... .......... 27% 34.6M 2s
    ##  37300K .......... .......... .......... .......... .......... 27% 39.9M 2s
    ##  37350K .......... .......... .......... .......... .......... 27%  144M 2s
    ##  37400K .......... .......... .......... .......... .......... 27% 98.9M 2s
    ##  37450K .......... .......... .......... .......... .......... 27%  156M 2s
    ##  37500K .......... .......... .......... .......... .......... 27%  139M 2s
    ##  37550K .......... .......... .......... .......... .......... 27%  169M 2s
    ##  37600K .......... .......... .......... .......... .......... 27%  151M 2s
    ##  37650K .......... .......... .......... .......... .......... 27%  104M 2s
    ##  37700K .......... .......... .......... .......... .......... 28%  105M 2s
    ##  37750K .......... .......... .......... .......... .......... 28%  103M 2s
    ##  37800K .......... .......... .......... .......... .......... 28%  102M 2s
    ##  37850K .......... .......... .......... .......... .......... 28%  146M 2s
    ##  37900K .......... .......... .......... .......... .......... 28%  109M 2s
    ##  37950K .......... .......... .......... .......... .......... 28% 93.4M 2s
    ##  38000K .......... .......... .......... .......... .......... 28% 89.5M 2s
    ##  38050K .......... .......... .......... .......... .......... 28%  133M 2s
    ##  38100K .......... .......... .......... .......... .......... 28%  142M 2s
    ##  38150K .......... .......... .......... .......... .......... 28%  160M 2s
    ##  38200K .......... .......... .......... .......... .......... 28% 53.3M 2s
    ##  38250K .......... .......... .......... .......... .......... 28%  148M 2s
    ##  38300K .......... .......... .......... .......... .......... 28%  108M 2s
    ##  38350K .......... .......... .......... .......... .......... 28% 62.0M 2s
    ##  38400K .......... .......... .......... .......... .......... 28% 80.7M 2s
    ##  38450K .......... .......... .......... .......... .......... 28%  159M 2s
    ##  38500K .......... .......... .......... .......... .......... 28%  140M 2s
    ##  38550K .......... .......... .......... .......... .......... 28% 31.7M 2s
    ##  38600K .......... .......... .......... .......... .......... 28% 76.7M 2s
    ##  38650K .......... .......... .......... .......... .......... 28% 27.4M 2s
    ##  38700K .......... .......... .......... .......... .......... 28% 50.2M 2s
    ##  38750K .......... .......... .......... .......... .......... 28% 47.7M 2s
    ##  38800K .......... .......... .......... .......... .......... 28%  151M 2s
    ##  38850K .......... .......... .......... .......... .......... 28% 44.4M 2s
    ##  38900K .......... .......... .......... .......... .......... 28%  156M 2s
    ##  38950K .......... .......... .......... .......... .......... 28%  139M 2s
    ##  39000K .......... .......... .......... .......... .......... 28% 42.3M 2s
    ##  39050K .......... .......... .......... .......... .......... 29%  116M 2s
    ##  39100K .......... .......... .......... .......... .......... 29% 62.2M 2s
    ##  39150K .......... .......... .......... .......... .......... 29% 84.3M 2s
    ##  39200K .......... .......... .......... .......... .......... 29%  102M 2s
    ##  39250K .......... .......... .......... .......... .......... 29%  104M 2s
    ##  39300K .......... .......... .......... .......... .......... 29%  162M 2s
    ##  39350K .......... .......... .......... .......... .......... 29% 27.0M 2s
    ##  39400K .......... .......... .......... .......... .......... 29%  103M 2s
    ##  39450K .......... .......... .......... .......... .......... 29% 53.3M 2s
    ##  39500K .......... .......... .......... .......... .......... 29% 87.9M 2s
    ##  39550K .......... .......... .......... .......... .......... 29% 69.3M 2s
    ##  39600K .......... .......... .......... .......... .......... 29% 57.7M 2s
    ##  39650K .......... .......... .......... .......... .......... 29% 83.2M 2s
    ##  39700K .......... .......... .......... .......... .......... 29%  104M 2s
    ##  39750K .......... .......... .......... .......... .......... 29%  119M 2s
    ##  39800K .......... .......... .......... .......... .......... 29%  108M 2s
    ##  39850K .......... .......... .......... .......... .......... 29% 52.7M 2s
    ##  39900K .......... .......... .......... .......... .......... 29% 57.2M 2s
    ##  39950K .......... .......... .......... .......... .......... 29% 58.5M 2s
    ##  40000K .......... .......... .......... .......... .......... 29% 81.4M 2s
    ##  40050K .......... .......... .......... .......... .......... 29%  113M 2s
    ##  40100K .......... .......... .......... .......... .......... 29% 68.9M 2s
    ##  40150K .......... .......... .......... .......... .......... 29% 79.5M 2s
    ##  40200K .......... .......... .......... .......... .......... 29% 83.4M 2s
    ##  40250K .......... .......... .......... .......... .......... 29%  102M 2s
    ##  40300K .......... .......... .......... .......... .......... 29% 56.0M 2s
    ##  40350K .......... .......... .......... .......... .......... 29% 50.3M 2s
    ##  40400K .......... .......... .......... .......... .......... 30% 77.1M 2s
    ##  40450K .......... .......... .......... .......... .......... 30% 46.0M 2s
    ##  40500K .......... .......... .......... .......... .......... 30%  124M 2s
    ##  40550K .......... .......... .......... .......... .......... 30%  101M 2s
    ##  40600K .......... .......... .......... .......... .......... 30% 75.8M 2s
    ##  40650K .......... .......... .......... .......... .......... 30% 72.4M 2s
    ##  40700K .......... .......... .......... .......... .......... 30% 81.3M 2s
    ##  40750K .......... .......... .......... .......... .......... 30% 90.6M 2s
    ##  40800K .......... .......... .......... .......... .......... 30% 76.1M 2s
    ##  40850K .......... .......... .......... .......... .......... 30% 71.8M 2s
    ##  40900K .......... .......... .......... .......... .......... 30% 70.4M 2s
    ##  40950K .......... .......... .......... .......... .......... 30%  121M 2s
    ##  41000K .......... .......... .......... .......... .......... 30% 61.8M 2s
    ##  41050K .......... .......... .......... .......... .......... 30% 75.2M 2s
    ##  41100K .......... .......... .......... .......... .......... 30% 62.9M 2s
    ##  41150K .......... .......... .......... .......... .......... 30%  104M 2s
    ##  41200K .......... .......... .......... .......... .......... 30% 73.9M 2s
    ##  41250K .......... .......... .......... .......... .......... 30% 86.0M 2s
    ##  41300K .......... .......... .......... .......... .......... 30% 64.1M 2s
    ##  41350K .......... .......... .......... .......... .......... 30% 59.2M 2s
    ##  41400K .......... .......... .......... .......... .......... 30%  118M 2s
    ##  41450K .......... .......... .......... .......... .......... 30%  107M 2s
    ##  41500K .......... .......... .......... .......... .......... 30% 70.9M 2s
    ##  41550K .......... .......... .......... .......... .......... 30% 67.6M 2s
    ##  41600K .......... .......... .......... .......... .......... 30% 86.9M 2s
    ##  41650K .......... .......... .......... .......... .......... 30%  110M 2s
    ##  41700K .......... .......... .......... .......... .......... 30% 64.8M 2s
    ##  41750K .......... .......... .......... .......... .......... 31% 73.5M 2s
    ##  41800K .......... .......... .......... .......... .......... 31% 79.3M 2s
    ##  41850K .......... .......... .......... .......... .......... 31% 58.9M 2s
    ##  41900K .......... .......... .......... .......... .......... 31% 73.4M 2s
    ##  41950K .......... .......... .......... .......... .......... 31%  167M 2s
    ##  42000K .......... .......... .......... .......... .......... 31% 81.8M 2s
    ##  42050K .......... .......... .......... .......... .......... 31% 87.2M 2s
    ##  42100K .......... .......... .......... .......... .......... 31% 69.1M 2s
    ##  42150K .......... .......... .......... .......... .......... 31%  155M 2s
    ##  42200K .......... .......... .......... .......... .......... 31% 73.5M 2s
    ##  42250K .......... .......... .......... .......... .......... 31% 65.8M 2s
    ##  42300K .......... .......... .......... .......... .......... 31% 92.9M 2s
    ##  42350K .......... .......... .......... .......... .......... 31%  137M 2s
    ##  42400K .......... .......... .......... .......... .......... 31% 58.2M 2s
    ##  42450K .......... .......... .......... .......... .......... 31% 91.4M 2s
    ##  42500K .......... .......... .......... .......... .......... 31% 86.3M 2s
    ##  42550K .......... .......... .......... .......... .......... 31%  116M 2s
    ##  42600K .......... .......... .......... .......... .......... 31% 72.5M 2s
    ##  42650K .......... .......... .......... .......... .......... 31% 97.5M 2s
    ##  42700K .......... .......... .......... .......... .......... 31% 91.7M 2s
    ##  42750K .......... .......... .......... .......... .......... 31% 79.3M 2s
    ##  42800K .......... .......... .......... .......... .......... 31%  102M 2s
    ##  42850K .......... .......... .......... .......... .......... 31%  101M 2s
    ##  42900K .......... .......... .......... .......... .......... 31%  105M 2s
    ##  42950K .......... .......... .......... .......... .......... 31% 75.9M 2s
    ##  43000K .......... .......... .......... .......... .......... 31% 84.1M 2s
    ##  43050K .......... .......... .......... .......... .......... 31%  107M 2s
    ##  43100K .......... .......... .......... .......... .......... 32% 77.9M 2s
    ##  43150K .......... .......... .......... .......... .......... 32%  103M 2s
    ##  43200K .......... .......... .......... .......... .......... 32% 86.5M 2s
    ##  43250K .......... .......... .......... .......... .......... 32%  131M 2s
    ##  43300K .......... .......... .......... .......... .......... 32% 72.7M 2s
    ##  43350K .......... .......... .......... .......... .......... 32%  101M 2s
    ##  43400K .......... .......... .......... .......... .......... 32% 76.5M 2s
    ##  43450K .......... .......... .......... .......... .......... 32%  116M 2s
    ##  43500K .......... .......... .......... .......... .......... 32% 71.4M 2s
    ##  43550K .......... .......... .......... .......... .......... 32%  182M 2s
    ##  43600K .......... .......... .......... .......... .......... 32% 88.4M 2s
    ##  43650K .......... .......... .......... .......... .......... 32% 67.2M 2s
    ##  43700K .......... .......... .......... .......... .......... 32%  106M 2s
    ##  43750K .......... .......... .......... .......... .......... 32% 69.5M 2s
    ##  43800K .......... .......... .......... .......... .......... 32%  102M 2s
    ##  43850K .......... .......... .......... .......... .......... 32%  107M 2s
    ##  43900K .......... .......... .......... .......... .......... 32% 91.1M 2s
    ##  43950K .......... .......... .......... .......... .......... 32% 79.7M 2s
    ##  44000K .......... .......... .......... .......... .......... 32% 75.1M 2s
    ##  44050K .......... .......... .......... .......... .......... 32% 80.7M 2s
    ##  44100K .......... .......... .......... .......... .......... 32% 85.0M 2s
    ##  44150K .......... .......... .......... .......... .......... 32%  189M 2s
    ##  44200K .......... .......... .......... .......... .......... 32% 69.0M 2s
    ##  44250K .......... .......... .......... .......... .......... 32% 96.9M 2s
    ##  44300K .......... .......... .......... .......... .......... 32% 88.3M 2s
    ##  44350K .......... .......... .......... .......... .......... 32% 90.0M 2s
    ##  44400K .......... .......... .......... .......... .......... 32% 97.8M 2s
    ##  44450K .......... .......... .......... .......... .......... 33% 96.3M 2s
    ##  44500K .......... .......... .......... .......... .......... 33% 90.9M 2s
    ##  44550K .......... .......... .......... .......... .......... 33% 89.8M 2s
    ##  44600K .......... .......... .......... .......... .......... 33% 86.8M 2s
    ##  44650K .......... .......... .......... .......... .......... 33% 99.5M 2s
    ##  44700K .......... .......... .......... .......... .......... 33% 73.6M 2s
    ##  44750K .......... .......... .......... .......... .......... 33% 94.5M 2s
    ##  44800K .......... .......... .......... .......... .......... 33% 71.4M 2s
    ##  44850K .......... .......... .......... .......... .......... 33% 79.5M 2s
    ##  44900K .......... .......... .......... .......... .......... 33% 73.4M 2s
    ##  44950K .......... .......... .......... .......... .......... 33%  122M 2s
    ##  45000K .......... .......... .......... .......... .......... 33% 45.0M 2s
    ##  45050K .......... .......... .......... .......... .......... 33% 95.6M 2s
    ##  45100K .......... .......... .......... .......... .......... 33% 81.8M 2s
    ##  45150K .......... .......... .......... .......... .......... 33% 88.3M 2s
    ##  45200K .......... .......... .......... .......... .......... 33% 80.3M 2s
    ##  45250K .......... .......... .......... .......... .......... 33% 97.7M 2s
    ##  45300K .......... .......... .......... .......... .......... 33% 64.4M 2s
    ##  45350K .......... .......... .......... .......... .......... 33%  104M 2s
    ##  45400K .......... .......... .......... .......... .......... 33% 98.9M 2s
    ##  45450K .......... .......... .......... .......... .......... 33%  128M 2s
    ##  45500K .......... .......... .......... .......... .......... 33%  101M 2s
    ##  45550K .......... .......... .......... .......... .......... 33%  111M 2s
    ##  45600K .......... .......... .......... .......... .......... 33%  145M 2s
    ##  45650K .......... .......... .......... .......... .......... 33%  123M 2s
    ##  45700K .......... .......... .......... .......... .......... 33% 77.1M 2s
    ##  45750K .......... .......... .......... .......... .......... 33%  129M 1s
    ##  45800K .......... .......... .......... .......... .......... 34%  122M 1s
    ##  45850K .......... .......... .......... .......... .......... 34%  146M 1s
    ##  45900K .......... .......... .......... .......... .......... 34%  148M 1s
    ##  45950K .......... .......... .......... .......... .......... 34%  164M 1s
    ##  46000K .......... .......... .......... .......... .......... 34%  106M 1s
    ##  46050K .......... .......... .......... .......... .......... 34%  143M 1s
    ##  46100K .......... .......... .......... .......... .......... 34%  149M 1s
    ##  46150K .......... .......... .......... .......... .......... 34%  116M 1s
    ##  46200K .......... .......... .......... .......... .......... 34%  130M 1s
    ##  46250K .......... .......... .......... .......... .......... 34%  144M 1s
    ##  46300K .......... .......... .......... .......... .......... 34% 97.2M 1s
    ##  46350K .......... .......... .......... .......... .......... 34%  150M 1s
    ##  46400K .......... .......... .......... .......... .......... 34%  158M 1s
    ##  46450K .......... .......... .......... .......... .......... 34% 55.0M 1s
    ##  46500K .......... .......... .......... .......... .......... 34%  124M 1s
    ##  46550K .......... .......... .......... .......... .......... 34%  137M 1s
    ##  46600K .......... .......... .......... .......... .......... 34%  158M 1s
    ##  46650K .......... .......... .......... .......... .......... 34% 94.4M 1s
    ##  46700K .......... .......... .......... .......... .......... 34% 94.4M 1s
    ##  46750K .......... .......... .......... .......... .......... 34% 81.9M 1s
    ##  46800K .......... .......... .......... .......... .......... 34%  147M 1s
    ##  46850K .......... .......... .......... .......... .......... 34% 79.4M 1s
    ##  46900K .......... .......... .......... .......... .......... 34%  111M 1s
    ##  46950K .......... .......... .......... .......... .......... 34% 76.4M 1s
    ##  47000K .......... .......... .......... .......... .......... 34% 48.3M 1s
    ##  47050K .......... .......... .......... .......... .......... 34% 64.0M 1s
    ##  47100K .......... .......... .......... .......... .......... 34%  104M 1s
    ##  47150K .......... .......... .......... .......... .......... 35% 25.8M 1s
    ##  47200K .......... .......... .......... .......... .......... 35%  146M 1s
    ##  47250K .......... .......... .......... .......... .......... 35%  156M 1s
    ##  47300K .......... .......... .......... .......... .......... 35%  147M 1s
    ##  47350K .......... .......... .......... .......... .......... 35%  172M 1s
    ##  47400K .......... .......... .......... .......... .......... 35% 61.9M 1s
    ##  47450K .......... .......... .......... .......... .......... 35%  140M 1s
    ##  47500K .......... .......... .......... .......... .......... 35% 86.2M 1s
    ##  47550K .......... .......... .......... .......... .......... 35%  125M 1s
    ##  47600K .......... .......... .......... .......... .......... 35% 68.6M 1s
    ##  47650K .......... .......... .......... .......... .......... 35% 54.2M 1s
    ##  47700K .......... .......... .......... .......... .......... 35% 36.1M 1s
    ##  47750K .......... .......... .......... .......... .......... 35% 65.0M 1s
    ##  47800K .......... .......... .......... .......... .......... 35% 65.7M 1s
    ##  47850K .......... .......... .......... .......... .......... 35%  135M 1s
    ##  47900K .......... .......... .......... .......... .......... 35% 51.8M 1s
    ##  47950K .......... .......... .......... .......... .......... 35%  134M 1s
    ##  48000K .......... .......... .......... .......... .......... 35%  158M 1s
    ##  48050K .......... .......... .......... .......... .......... 35%  114M 1s
    ##  48100K .......... .......... .......... .......... .......... 35%  165M 1s
    ##  48150K .......... .......... .......... .......... .......... 35% 20.5M 1s
    ##  48200K .......... .......... .......... .......... .......... 35%  140M 1s
    ##  48250K .......... .......... .......... .......... .......... 35% 74.4M 1s
    ##  48300K .......... .......... .......... .......... .......... 35% 54.9M 1s
    ##  48350K .......... .......... .......... .......... .......... 35%  149M 1s
    ##  48400K .......... .......... .......... .......... .......... 35%  148M 1s
    ##  48450K .......... .......... .......... .......... .......... 35%  168M 1s
    ##  48500K .......... .......... .......... .......... .......... 36% 41.3M 1s
    ##  48550K .......... .......... .......... .......... .......... 36% 75.5M 1s
    ##  48600K .......... .......... .......... .......... .......... 36%  106M 1s
    ##  48650K .......... .......... .......... .......... .......... 36% 74.3M 1s
    ##  48700K .......... .......... .......... .......... .......... 36% 70.5M 1s
    ##  48750K .......... .......... .......... .......... .......... 36%  131M 1s
    ##  48800K .......... .......... .......... .......... .......... 36%  151M 1s
    ##  48850K .......... .......... .......... .......... .......... 36%  150M 1s
    ##  48900K .......... .......... .......... .......... .......... 36% 77.1M 1s
    ##  48950K .......... .......... .......... .......... .......... 36% 68.0M 1s
    ##  49000K .......... .......... .......... .......... .......... 36% 96.3M 1s
    ##  49050K .......... .......... .......... .......... .......... 36% 31.2M 1s
    ##  49100K .......... .......... .......... .......... .......... 36% 40.2M 1s
    ##  49150K .......... .......... .......... .......... .......... 36% 48.1M 1s
    ##  49200K .......... .......... .......... .......... .......... 36%  144M 1s
    ##  49250K .......... .......... .......... .......... .......... 36%  129M 1s
    ##  49300K .......... .......... .......... .......... .......... 36% 40.9M 1s
    ##  49350K .......... .......... .......... .......... .......... 36%  169M 1s
    ##  49400K .......... .......... .......... .......... .......... 36%  132M 1s
    ##  49450K .......... .......... .......... .......... .......... 36%  156M 1s
    ##  49500K .......... .......... .......... .......... .......... 36% 48.9M 1s
    ##  49550K .......... .......... .......... .......... .......... 36%  105M 1s
    ##  49600K .......... .......... .......... .......... .......... 36% 55.7M 1s
    ##  49650K .......... .......... .......... .......... .......... 36%  171M 1s
    ##  49700K .......... .......... .......... .......... .......... 36% 31.8M 1s
    ##  49750K .......... .......... .......... .......... .......... 36% 66.8M 1s
    ##  49800K .......... .......... .......... .......... .......... 36%  130M 1s
    ##  49850K .......... .......... .......... .......... .......... 37% 60.9M 1s
    ##  49900K .......... .......... .......... .......... .......... 37% 90.3M 1s
    ##  49950K .......... .......... .......... .......... .......... 37%  183M 1s
    ##  50000K .......... .......... .......... .......... .......... 37%  113M 1s
    ##  50050K .......... .......... .......... .......... .......... 37%  144M 1s
    ##  50100K .......... .......... .......... .......... .......... 37%  115M 1s
    ##  50150K .......... .......... .......... .......... .......... 37%  160M 1s
    ##  50200K .......... .......... .......... .......... .......... 37% 30.3M 1s
    ##  50250K .......... .......... .......... .......... .......... 37% 92.0M 1s
    ##  50300K .......... .......... .......... .......... .......... 37% 70.3M 1s
    ##  50350K .......... .......... .......... .......... .......... 37% 85.1M 1s
    ##  50400K .......... .......... .......... .......... .......... 37%  144M 1s
    ##  50450K .......... .......... .......... .......... .......... 37%  159M 1s
    ##  50500K .......... .......... .......... .......... .......... 37%  104M 1s
    ##  50550K .......... .......... .......... .......... .......... 37%  123M 1s
    ##  50600K .......... .......... .......... .......... .......... 37%  144M 1s
    ##  50650K .......... .......... .......... .......... .......... 37% 34.5M 1s
    ##  50700K .......... .......... .......... .......... .......... 37%  138M 1s
    ##  50750K .......... .......... .......... .......... .......... 37% 17.6M 1s
    ##  50800K .......... .......... .......... .......... .......... 37%  114M 1s
    ##  50850K .......... .......... .......... .......... .......... 37%  166M 1s
    ##  50900K .......... .......... .......... .......... .......... 37%  137M 1s
    ##  50950K .......... .......... .......... .......... .......... 37%  123M 1s
    ##  51000K .......... .......... .......... .......... .......... 37% 41.2M 1s
    ##  51050K .......... .......... .......... .......... .......... 37% 63.9M 1s
    ##  51100K .......... .......... .......... .......... .......... 37% 26.4M 1s
    ##  51150K .......... .......... .......... .......... .......... 37% 60.6M 1s
    ##  51200K .......... .......... .......... .......... .......... 38% 41.2M 1s
    ##  51250K .......... .......... .......... .......... .......... 38%  151M 1s
    ##  51300K .......... .......... .......... .......... .......... 38% 19.3M 1s
    ##  51350K .......... .......... .......... .......... .......... 38%  149M 1s
    ##  51400K .......... .......... .......... .......... .......... 38%  136M 1s
    ##  51450K .......... .......... .......... .......... .......... 38%  131M 1s
    ##  51500K .......... .......... .......... .......... .......... 38% 85.9M 1s
    ##  51550K .......... .......... .......... .......... .......... 38%  106M 1s
    ##  51600K .......... .......... .......... .......... .......... 38%  125M 1s
    ##  51650K .......... .......... .......... .......... .......... 38%  149M 1s
    ##  51700K .......... .......... .......... .......... .......... 38%  142M 1s
    ##  51750K .......... .......... .......... .......... .......... 38% 63.8M 1s
    ##  51800K .......... .......... .......... .......... .......... 38% 73.7M 1s
    ##  51850K .......... .......... .......... .......... .......... 38% 97.8M 1s
    ##  51900K .......... .......... .......... .......... .......... 38%  145M 1s
    ##  51950K .......... .......... .......... .......... .......... 38%  156M 1s
    ##  52000K .......... .......... .......... .......... .......... 38%  140M 1s
    ##  52050K .......... .......... .......... .......... .......... 38%  180M 1s
    ##  52100K .......... .......... .......... .......... .......... 38% 39.4M 1s
    ##  52150K .......... .......... .......... .......... .......... 38%  114M 1s
    ##  52200K .......... .......... .......... .......... .......... 38% 81.2M 1s
    ##  52250K .......... .......... .......... .......... .......... 38% 82.4M 1s
    ##  52300K .......... .......... .......... .......... .......... 38% 83.0M 1s
    ##  52350K .......... .......... .......... .......... .......... 38% 82.0M 1s
    ##  52400K .......... .......... .......... .......... .......... 38%  147M 1s
    ##  52450K .......... .......... .......... .......... .......... 38%  106M 1s
    ##  52500K .......... .......... .......... .......... .......... 39% 91.9M 1s
    ##  52550K .......... .......... .......... .......... .......... 39%  114M 1s
    ##  52600K .......... .......... .......... .......... .......... 39% 82.5M 1s
    ##  52650K .......... .......... .......... .......... .......... 39%  167M 1s
    ##  52700K .......... .......... .......... .......... .......... 39% 80.6M 1s
    ##  52750K .......... .......... .......... .......... .......... 39% 81.9M 1s
    ##  52800K .......... .......... .......... .......... .......... 39%  101M 1s
    ##  52850K .......... .......... .......... .......... .......... 39% 69.0M 1s
    ##  52900K .......... .......... .......... .......... .......... 39% 67.9M 1s
    ##  52950K .......... .......... .......... .......... .......... 39%  104M 1s
    ##  53000K .......... .......... .......... .......... .......... 39% 74.4M 1s
    ##  53050K .......... .......... .......... .......... .......... 39%  111M 1s
    ##  53100K .......... .......... .......... .......... .......... 39%  130M 1s
    ##  53150K .......... .......... .......... .......... .......... 39%  126M 1s
    ##  53200K .......... .......... .......... .......... .......... 39% 83.9M 1s
    ##  53250K .......... .......... .......... .......... .......... 39%  123M 1s
    ##  53300K .......... .......... .......... .......... .......... 39% 54.0M 1s
    ##  53350K .......... .......... .......... .......... .......... 39% 56.5M 1s
    ##  53400K .......... .......... .......... .......... .......... 39% 91.1M 1s
    ##  53450K .......... .......... .......... .......... .......... 39%  128M 1s
    ##  53500K .......... .......... .......... .......... .......... 39%  141M 1s
    ##  53550K .......... .......... .......... .......... .......... 39% 63.9M 1s
    ##  53600K .......... .......... .......... .......... .......... 39%  145M 1s
    ##  53650K .......... .......... .......... .......... .......... 39%  117M 1s
    ##  53700K .......... .......... .......... .......... .......... 39%  128M 1s
    ##  53750K .......... .......... .......... .......... .......... 39% 96.8M 1s
    ##  53800K .......... .......... .......... .......... .......... 39% 54.4M 1s
    ##  53850K .......... .......... .......... .......... .......... 40% 91.9M 1s
    ##  53900K .......... .......... .......... .......... .......... 40% 31.4M 1s
    ##  53950K .......... .......... .......... .......... .......... 40% 55.4M 1s
    ##  54000K .......... .......... .......... .......... .......... 40% 34.7M 1s
    ##  54050K .......... .......... .......... .......... .......... 40%  155M 1s
    ##  54100K .......... .......... .......... .......... .......... 40%  144M 1s
    ##  54150K .......... .......... .......... .......... .......... 40%  149M 1s
    ##  54200K .......... .......... .......... .......... .......... 40%  150M 1s
    ##  54250K .......... .......... .......... .......... .......... 40% 44.0M 1s
    ##  54300K .......... .......... .......... .......... .......... 40%  138M 1s
    ##  54350K .......... .......... .......... .......... .......... 40%  118M 1s
    ##  54400K .......... .......... .......... .......... .......... 40% 91.0M 1s
    ##  54450K .......... .......... .......... .......... .......... 40%  153M 1s
    ##  54500K .......... .......... .......... .......... .......... 40% 80.0M 1s
    ##  54550K .......... .......... .......... .......... .......... 40% 88.8M 1s
    ##  54600K .......... .......... .......... .......... .......... 40%  103M 1s
    ##  54650K .......... .......... .......... .......... .......... 40% 73.2M 1s
    ##  54700K .......... .......... .......... .......... .......... 40% 74.6M 1s
    ##  54750K .......... .......... .......... .......... .......... 40% 70.2M 1s
    ##  54800K .......... .......... .......... .......... .......... 40% 72.6M 1s
    ##  54850K .......... .......... .......... .......... .......... 40% 92.7M 1s
    ##  54900K .......... .......... .......... .......... .......... 40% 35.6M 1s
    ##  54950K .......... .......... .......... .......... .......... 40%  124M 1s
    ##  55000K .......... .......... .......... .......... .......... 40%  139M 1s
    ##  55050K .......... .......... .......... .......... .......... 40%  164M 1s
    ##  55100K .......... .......... .......... .......... .......... 40% 12.2M 1s
    ##  55150K .......... .......... .......... .......... .......... 40%  101M 1s
    ##  55200K .......... .......... .......... .......... .......... 41%  148M 1s
    ##  55250K .......... .......... .......... .......... .......... 41%  142M 1s
    ##  55300K .......... .......... .......... .......... .......... 41%  148M 1s
    ##  55350K .......... .......... .......... .......... .......... 41%  147M 1s
    ##  55400K .......... .......... .......... .......... .......... 41%  150M 1s
    ##  55450K .......... .......... .......... .......... .......... 41%  167M 1s
    ##  55500K .......... .......... .......... .......... .......... 41% 30.0M 1s
    ##  55550K .......... .......... .......... .......... .......... 41% 58.9M 1s
    ##  55600K .......... .......... .......... .......... .......... 41% 23.5M 1s
    ##  55650K .......... .......... .......... .......... .......... 41%  130M 1s
    ##  55700K .......... .......... .......... .......... .......... 41% 49.8M 1s
    ##  55750K .......... .......... .......... .......... .......... 41% 57.1M 1s
    ##  55800K .......... .......... .......... .......... .......... 41%  112M 1s
    ##  55850K .......... .......... .......... .......... .......... 41% 56.8M 1s
    ##  55900K .......... .......... .......... .......... .......... 41% 54.0M 1s
    ##  55950K .......... .......... .......... .......... .......... 41%  110M 1s
    ##  56000K .......... .......... .......... .......... .......... 41%  118M 1s
    ##  56050K .......... .......... .......... .......... .......... 41% 85.0M 1s
    ##  56100K .......... .......... .......... .......... .......... 41%  121M 1s
    ##  56150K .......... .......... .......... .......... .......... 41%  117M 1s
    ##  56200K .......... .......... .......... .......... .......... 41%  132M 1s
    ##  56250K .......... .......... .......... .......... .......... 41% 12.3M 1s
    ##  56300K .......... .......... .......... .......... .......... 41%  106M 1s
    ##  56350K .......... .......... .......... .......... .......... 41%  163M 1s
    ##  56400K .......... .......... .......... .......... .......... 41%  145M 1s
    ##  56450K .......... .......... .......... .......... .......... 41%  152M 1s
    ##  56500K .......... .......... .......... .......... .......... 41%  137M 1s
    ##  56550K .......... .......... .......... .......... .......... 42%  170M 1s
    ##  56600K .......... .......... .......... .......... .......... 42%  133M 1s
    ##  56650K .......... .......... .......... .......... .......... 42% 67.9M 1s
    ##  56700K .......... .......... .......... .......... .......... 42% 40.2M 1s
    ##  56750K .......... .......... .......... .......... .......... 42% 24.6M 1s
    ##  56800K .......... .......... .......... .......... .......... 42% 39.9M 1s
    ##  56850K .......... .......... .......... .......... .......... 42% 90.8M 1s
    ##  56900K .......... .......... .......... .......... .......... 42%  121M 1s
    ##  56950K .......... .......... .......... .......... .......... 42% 87.8M 1s
    ##  57000K .......... .......... .......... .......... .......... 42% 75.8M 1s
    ##  57050K .......... .......... .......... .......... .......... 42% 77.0M 1s
    ##  57100K .......... .......... .......... .......... .......... 42%  133M 1s
    ##  57150K .......... .......... .......... .......... .......... 42%  166M 1s
    ##  57200K .......... .......... .......... .......... .......... 42%  145M 1s
    ##  57250K .......... .......... .......... .......... .......... 42%  171M 1s
    ##  57300K .......... .......... .......... .......... .......... 42% 63.3M 1s
    ##  57350K .......... .......... .......... .......... .......... 42% 56.7M 1s
    ##  57400K .......... .......... .......... .......... .......... 42% 65.6M 1s
    ##  57450K .......... .......... .......... .......... .......... 42%  155M 1s
    ##  57500K .......... .......... .......... .......... .......... 42%  132M 1s
    ##  57550K .......... .......... .......... .......... .......... 42%  179M 1s
    ##  57600K .......... .......... .......... .......... .......... 42%  136M 1s
    ##  57650K .......... .......... .......... .......... .......... 42% 74.4M 1s
    ##  57700K .......... .......... .......... .......... .......... 42% 60.4M 1s
    ##  57750K .......... .......... .......... .......... .......... 42%  149M 1s
    ##  57800K .......... .......... .......... .......... .......... 42% 16.2M 1s
    ##  57850K .......... .......... .......... .......... .......... 42% 21.9M 1s
    ##  57900K .......... .......... .......... .......... .......... 43% 56.3M 1s
    ##  57950K .......... .......... .......... .......... .......... 43% 48.9M 1s
    ##  58000K .......... .......... .......... .......... .......... 43% 41.7M 1s
    ##  58050K .......... .......... .......... .......... .......... 43% 72.0M 1s
    ##  58100K .......... .......... .......... .......... .......... 43% 66.3M 1s
    ##  58150K .......... .......... .......... .......... .......... 43% 89.4M 1s
    ##  58200K .......... .......... .......... .......... .......... 43% 68.2M 1s
    ##  58250K .......... .......... .......... .......... .......... 43%  151M 1s
    ##  58300K .......... .......... .......... .......... .......... 43% 83.6M 1s
    ##  58350K .......... .......... .......... .......... .......... 43%  167M 1s
    ##  58400K .......... .......... .......... .......... .......... 43%  140M 1s
    ##  58450K .......... .......... .......... .......... .......... 43% 13.5M 1s
    ##  58500K .......... .......... .......... .......... .......... 43%  131M 1s
    ##  58550K .......... .......... .......... .......... .......... 43%  147M 1s
    ##  58600K .......... .......... .......... .......... .......... 43%  140M 1s
    ##  58650K .......... .......... .......... .......... .......... 43%  139M 1s
    ##  58700K .......... .......... .......... .......... .......... 43%  120M 1s
    ##  58750K .......... .......... .......... .......... .......... 43% 88.8M 1s
    ##  58800K .......... .......... .......... .......... .......... 43% 88.5M 1s
    ##  58850K .......... .......... .......... .......... .......... 43% 72.2M 1s
    ##  58900K .......... .......... .......... .......... .......... 43% 32.9M 1s
    ##  58950K .......... .......... .......... .......... .......... 43% 38.9M 1s
    ##  59000K .......... .......... .......... .......... .......... 43% 94.5M 1s
    ##  59050K .......... .......... .......... .......... .......... 43%  178M 1s
    ##  59100K .......... .......... .......... .......... .......... 43% 85.1M 1s
    ##  59150K .......... .......... .......... .......... .......... 43% 96.5M 1s
    ##  59200K .......... .......... .......... .......... .......... 43% 82.7M 1s
    ##  59250K .......... .......... .......... .......... .......... 44% 82.1M 1s
    ##  59300K .......... .......... .......... .......... .......... 44% 81.4M 1s
    ##  59350K .......... .......... .......... .......... .......... 44%  144M 1s
    ##  59400K .......... .......... .......... .......... .......... 44%  134M 1s
    ##  59450K .......... .......... .......... .......... .......... 44% 15.1M 1s
    ##  59500K .......... .......... .......... .......... .......... 44%  136M 1s
    ##  59550K .......... .......... .......... .......... .......... 44% 87.8M 1s
    ##  59600K .......... .......... .......... .......... .......... 44%  158M 1s
    ##  59650K .......... .......... .......... .......... .......... 44%  129M 1s
    ##  59700K .......... .......... .......... .......... .......... 44%  127M 1s
    ##  59750K .......... .......... .......... .......... .......... 44%  123M 1s
    ##  59800K .......... .......... .......... .......... .......... 44%  111M 1s
    ##  59850K .......... .......... .......... .......... .......... 44%  105M 1s
    ##  59900K .......... .......... .......... .......... .......... 44%  109M 1s
    ##  59950K .......... .......... .......... .......... .......... 44%  137M 1s
    ##  60000K .......... .......... .......... .......... .......... 44% 31.8M 1s
    ##  60050K .......... .......... .......... .......... .......... 44% 83.0M 1s
    ##  60100K .......... .......... .......... .......... .......... 44% 64.4M 1s
    ##  60150K .......... .......... .......... .......... .......... 44%  117M 1s
    ##  60200K .......... .......... .......... .......... .......... 44% 77.9M 1s
    ##  60250K .......... .......... .......... .......... .......... 44%  101M 1s
    ##  60300K .......... .......... .......... .......... .......... 44% 82.2M 1s
    ##  60350K .......... .......... .......... .......... .......... 44% 84.7M 1s
    ##  60400K .......... .......... .......... .......... .......... 44%  121M 1s
    ##  60450K .......... .......... .......... .......... .......... 44% 66.6M 1s
    ##  60500K .......... .......... .......... .......... .......... 44%  105M 1s
    ##  60550K .......... .......... .......... .......... .......... 44% 97.6M 1s
    ##  60600K .......... .......... .......... .......... .......... 45% 43.8M 1s
    ##  60650K .......... .......... .......... .......... .......... 45% 83.2M 1s
    ##  60700K .......... .......... .......... .......... .......... 45% 85.7M 1s
    ##  60750K .......... .......... .......... .......... .......... 45% 84.6M 1s
    ##  60800K .......... .......... .......... .......... .......... 45% 99.9M 1s
    ##  60850K .......... .......... .......... .......... .......... 45% 49.3M 1s
    ##  60900K .......... .......... .......... .......... .......... 45% 85.1M 1s
    ##  60950K .......... .......... .......... .......... .......... 45% 87.7M 1s
    ##  61000K .......... .......... .......... .......... .......... 45% 74.0M 1s
    ##  61050K .......... .......... .......... .......... .......... 45%  133M 1s
    ##  61100K .......... .......... .......... .......... .......... 45% 74.0M 1s
    ##  61150K .......... .......... .......... .......... .......... 45% 83.8M 1s
    ##  61200K .......... .......... .......... .......... .......... 45%  120M 1s
    ##  61250K .......... .......... .......... .......... .......... 45% 48.3M 1s
    ##  61300K .......... .......... .......... .......... .......... 45%  138M 1s
    ##  61350K .......... .......... .......... .......... .......... 45% 95.8M 1s
    ##  61400K .......... .......... .......... .......... .......... 45% 71.2M 1s
    ##  61450K .......... .......... .......... .......... .......... 45%  101M 1s
    ##  61500K .......... .......... .......... .......... .......... 45% 92.5M 1s
    ##  61550K .......... .......... .......... .......... .......... 45% 74.8M 1s
    ##  61600K .......... .......... .......... .......... .......... 45% 81.1M 1s
    ##  61650K .......... .......... .......... .......... .......... 45% 63.5M 1s
    ##  61700K .......... .......... .......... .......... .......... 45%  106M 1s
    ##  61750K .......... .......... .......... .......... .......... 45%  106M 1s
    ##  61800K .......... .......... .......... .......... .......... 45% 74.7M 1s
    ##  61850K .......... .......... .......... .......... .......... 45% 97.7M 1s
    ##  61900K .......... .......... .......... .......... .......... 45% 60.3M 1s
    ##  61950K .......... .......... .......... .......... .......... 46%  147M 1s
    ##  62000K .......... .......... .......... .......... .......... 46% 72.5M 1s
    ##  62050K .......... .......... .......... .......... .......... 46% 77.3M 1s
    ##  62100K .......... .......... .......... .......... .......... 46% 73.8M 1s
    ##  62150K .......... .......... .......... .......... .......... 46% 79.8M 1s
    ##  62200K .......... .......... .......... .......... .......... 46%  122M 1s
    ##  62250K .......... .......... .......... .......... .......... 46% 77.2M 1s
    ##  62300K .......... .......... .......... .......... .......... 46% 89.1M 1s
    ##  62350K .......... .......... .......... .......... .......... 46% 83.4M 1s
    ##  62400K .......... .......... .......... .......... .......... 46% 83.9M 1s
    ##  62450K .......... .......... .......... .......... .......... 46% 77.0M 1s
    ##  62500K .......... .......... .......... .......... .......... 46% 76.0M 1s
    ##  62550K .......... .......... .......... .......... .......... 46%  143M 1s
    ##  62600K .......... .......... .......... .......... .......... 46% 11.1M 1s
    ##  62650K .......... .......... .......... .......... .......... 46% 39.0M 1s
    ##  62700K .......... .......... .......... .......... .......... 46% 39.4M 1s
    ##  62750K .......... .......... .......... .......... .......... 46% 91.8M 1s
    ##  62800K .......... .......... .......... .......... .......... 46% 36.5M 1s
    ##  62850K .......... .......... .......... .......... .......... 46%  146M 1s
    ##  62900K .......... .......... .......... .......... .......... 46%  150M 1s
    ##  62950K .......... .......... .......... .......... .......... 46%  154M 1s
    ##  63000K .......... .......... .......... .......... .......... 46%  118M 1s
    ##  63050K .......... .......... .......... .......... .......... 46%  124M 1s
    ##  63100K .......... .......... .......... .......... .......... 46% 81.1M 1s
    ##  63150K .......... .......... .......... .......... .......... 46% 90.2M 1s
    ##  63200K .......... .......... .......... .......... .......... 46% 81.2M 1s
    ##  63250K .......... .......... .......... .......... .......... 46%  113M 1s
    ##  63300K .......... .......... .......... .......... .......... 47% 26.9M 1s
    ##  63350K .......... .......... .......... .......... .......... 47% 81.7M 1s
    ##  63400K .......... .......... .......... .......... .......... 47% 75.6M 1s
    ##  63450K .......... .......... .......... .......... .......... 47% 64.2M 1s
    ##  63500K .......... .......... .......... .......... .......... 47%  146M 1s
    ##  63550K .......... .......... .......... .......... .......... 47%  117M 1s
    ##  63600K .......... .......... .......... .......... .......... 47% 95.3M 1s
    ##  63650K .......... .......... .......... .......... .......... 47%  129M 1s
    ##  63700K .......... .......... .......... .......... .......... 47%  119M 1s
    ##  63750K .......... .......... .......... .......... .......... 47%  112M 1s
    ##  63800K .......... .......... .......... .......... .......... 47% 83.9M 1s
    ##  63850K .......... .......... .......... .......... .......... 47% 61.2M 1s
    ##  63900K .......... .......... .......... .......... .......... 47%  102M 1s
    ##  63950K .......... .......... .......... .......... .......... 47% 88.5M 1s
    ##  64000K .......... .......... .......... .......... .......... 47% 41.4M 1s
    ##  64050K .......... .......... .......... .......... .......... 47% 55.5M 1s
    ##  64100K .......... .......... .......... .......... .......... 47% 87.6M 1s
    ##  64150K .......... .......... .......... .......... .......... 47%  133M 1s
    ##  64200K .......... .......... .......... .......... .......... 47% 47.3M 1s
    ##  64250K .......... .......... .......... .......... .......... 47% 64.2M 1s
    ##  64300K .......... .......... .......... .......... .......... 47% 89.0M 1s
    ##  64350K .......... .......... .......... .......... .......... 47%  136M 1s
    ##  64400K .......... .......... .......... .......... .......... 47%  103M 1s
    ##  64450K .......... .......... .......... .......... .......... 47% 60.7M 1s
    ##  64500K .......... .......... .......... .......... .......... 47% 65.1M 1s
    ##  64550K .......... .......... .......... .......... .......... 47%  112M 1s
    ##  64600K .......... .......... .......... .......... .......... 47% 95.3M 1s
    ##  64650K .......... .......... .......... .......... .......... 48% 44.2M 1s
    ##  64700K .......... .......... .......... .......... .......... 48%  149M 1s
    ##  64750K .......... .......... .......... .......... .......... 48% 84.4M 1s
    ##  64800K .......... .......... .......... .......... .......... 48% 78.6M 1s
    ##  64850K .......... .......... .......... .......... .......... 48% 86.3M 1s
    ##  64900K .......... .......... .......... .......... .......... 48%  111M 1s
    ##  64950K .......... .......... .......... .......... .......... 48% 75.0M 1s
    ##  65000K .......... .......... .......... .......... .......... 48% 97.7M 1s
    ##  65050K .......... .......... .......... .......... .......... 48% 52.7M 1s
    ##  65100K .......... .......... .......... .......... .......... 48% 59.1M 1s
    ##  65150K .......... .......... .......... .......... .......... 48% 59.1M 1s
    ##  65200K .......... .......... .......... .......... .......... 48% 91.5M 1s
    ##  65250K .......... .......... .......... .......... .......... 48%  153M 1s
    ##  65300K .......... .......... .......... .......... .......... 48% 43.9M 1s
    ##  65350K .......... .......... .......... .......... .......... 48%  153M 1s
    ##  65400K .......... .......... .......... .......... .......... 48% 68.1M 1s
    ##  65450K .......... .......... .......... .......... .......... 48% 35.7M 1s
    ##  65500K .......... .......... .......... .......... .......... 48% 33.9M 1s
    ##  65550K .......... .......... .......... .......... .......... 48% 37.9M 1s
    ##  65600K .......... .......... .......... .......... .......... 48% 42.9M 1s
    ##  65650K .......... .......... .......... .......... .......... 48% 37.5M 1s
    ##  65700K .......... .......... .......... .......... .......... 48% 41.6M 1s
    ##  65750K .......... .......... .......... .......... .......... 48% 25.1M 1s
    ##  65800K .......... .......... .......... .......... .......... 48% 39.7M 1s
    ##  65850K .......... .......... .......... .......... .......... 48% 94.0M 1s
    ##  65900K .......... .......... .......... .......... .......... 48% 84.8M 1s
    ##  65950K .......... .......... .......... .......... .......... 48% 97.6M 1s
    ##  66000K .......... .......... .......... .......... .......... 49% 66.2M 1s
    ##  66050K .......... .......... .......... .......... .......... 49%  126M 1s
    ##  66100K .......... .......... .......... .......... .......... 49%  100M 1s
    ##  66150K .......... .......... .......... .......... .......... 49% 80.9M 1s
    ##  66200K .......... .......... .......... .......... .......... 49% 99.8M 1s
    ##  66250K .......... .......... .......... .......... .......... 49%  105M 1s
    ##  66300K .......... .......... .......... .......... .......... 49% 82.8M 1s
    ##  66350K .......... .......... .......... .......... .......... 49% 75.6M 1s
    ##  66400K .......... .......... .......... .......... .......... 49% 28.8M 1s
    ##  66450K .......... .......... .......... .......... .......... 49%  103M 1s
    ##  66500K .......... .......... .......... .......... .......... 49% 92.2M 1s
    ##  66550K .......... .......... .......... .......... .......... 49% 66.6M 1s
    ##  66600K .......... .......... .......... .......... .......... 49% 74.6M 1s
    ##  66650K .......... .......... .......... .......... .......... 49%  102M 1s
    ##  66700K .......... .......... .......... .......... .......... 49%  116M 1s
    ##  66750K .......... .......... .......... .......... .......... 49% 91.6M 1s
    ##  66800K .......... .......... .......... .......... .......... 49% 75.8M 1s
    ##  66850K .......... .......... .......... .......... .......... 49%  101M 1s
    ##  66900K .......... .......... .......... .......... .......... 49% 65.5M 1s
    ##  66950K .......... .......... .......... .......... .......... 49% 80.1M 1s
    ##  67000K .......... .......... .......... .......... .......... 49% 75.9M 1s
    ##  67050K .......... .......... .......... .......... .......... 49%  107M 1s
    ##  67100K .......... .......... .......... .......... .......... 49% 94.4M 1s
    ##  67150K .......... .......... .......... .......... .......... 49% 85.2M 1s
    ##  67200K .......... .......... .......... .......... .......... 49% 76.1M 1s
    ##  67250K .......... .......... .......... .......... .......... 49% 68.8M 1s
    ##  67300K .......... .......... .......... .......... .......... 49% 21.3M 1s
    ##  67350K .......... .......... .......... .......... .......... 50% 56.2M 1s
    ##  67400K .......... .......... .......... .......... .......... 50% 65.2M 1s
    ##  67450K .......... .......... .......... .......... .......... 50% 74.9M 1s
    ##  67500K .......... .......... .......... .......... .......... 50% 89.2M 1s
    ##  67550K .......... .......... .......... .......... .......... 50% 73.5M 1s
    ##  67600K .......... .......... .......... .......... .......... 50% 88.4M 1s
    ##  67650K .......... .......... .......... .......... .......... 50% 97.9M 1s
    ##  67700K .......... .......... .......... .......... .......... 50% 82.0M 1s
    ##  67750K .......... .......... .......... .......... .......... 50% 83.5M 1s
    ##  67800K .......... .......... .......... .......... .......... 50% 52.5M 1s
    ##  67850K .......... .......... .......... .......... .......... 50% 73.2M 1s
    ##  67900K .......... .......... .......... .......... .......... 50% 76.2M 1s
    ##  67950K .......... .......... .......... .......... .......... 50% 83.3M 1s
    ##  68000K .......... .......... .......... .......... .......... 50% 69.4M 1s
    ##  68050K .......... .......... .......... .......... .......... 50% 63.4M 1s
    ##  68100K .......... .......... .......... .......... .......... 50% 70.3M 1s
    ##  68150K .......... .......... .......... .......... .......... 50% 61.6M 1s
    ##  68200K .......... .......... .......... .......... .......... 50% 66.2M 1s
    ##  68250K .......... .......... .......... .......... .......... 50% 92.4M 1s
    ##  68300K .......... .......... .......... .......... .......... 50%  100M 1s
    ##  68350K .......... .......... .......... .......... .......... 50% 88.7M 1s
    ##  68400K .......... .......... .......... .......... .......... 50% 78.9M 1s
    ##  68450K .......... .......... .......... .......... .......... 50% 73.3M 1s
    ##  68500K .......... .......... .......... .......... .......... 50% 68.5M 1s
    ##  68550K .......... .......... .......... .......... .......... 50% 88.4M 1s
    ##  68600K .......... .......... .......... .......... .......... 50% 69.6M 1s
    ##  68650K .......... .......... .......... .......... .......... 50% 74.5M 1s
    ##  68700K .......... .......... .......... .......... .......... 51% 79.0M 1s
    ##  68750K .......... .......... .......... .......... .......... 51% 94.7M 1s
    ##  68800K .......... .......... .......... .......... .......... 51% 38.6M 1s
    ##  68850K .......... .......... .......... .......... .......... 51%  110M 1s
    ##  68900K .......... .......... .......... .......... .......... 51% 95.1M 1s
    ##  68950K .......... .......... .......... .......... .......... 51% 87.5M 1s
    ##  69000K .......... .......... .......... .......... .......... 51% 83.4M 1s
    ##  69050K .......... .......... .......... .......... .......... 51% 35.5M 1s
    ##  69100K .......... .......... .......... .......... .......... 51% 59.7M 1s
    ##  69150K .......... .......... .......... .......... .......... 51%  104M 1s
    ##  69200K .......... .......... .......... .......... .......... 51% 94.7M 1s
    ##  69250K .......... .......... .......... .......... .......... 51% 96.6M 1s
    ##  69300K .......... .......... .......... .......... .......... 51% 83.0M 1s
    ##  69350K .......... .......... .......... .......... .......... 51% 78.3M 1s
    ##  69400K .......... .......... .......... .......... .......... 51% 77.1M 1s
    ##  69450K .......... .......... .......... .......... .......... 51% 94.3M 1s
    ##  69500K .......... .......... .......... .......... .......... 51% 72.0M 1s
    ##  69550K .......... .......... .......... .......... .......... 51%  101M 1s
    ##  69600K .......... .......... .......... .......... .......... 51% 85.9M 1s
    ##  69650K .......... .......... .......... .......... .......... 51% 93.3M 1s
    ##  69700K .......... .......... .......... .......... .......... 51% 92.2M 1s
    ##  69750K .......... .......... .......... .......... .......... 51%  109M 1s
    ##  69800K .......... .......... .......... .......... .......... 51% 94.4M 1s
    ##  69850K .......... .......... .......... .......... .......... 51% 83.2M 1s
    ##  69900K .......... .......... .......... .......... .......... 51% 82.9M 1s
    ##  69950K .......... .......... .......... .......... .......... 51% 90.7M 1s
    ##  70000K .......... .......... .......... .......... .......... 51% 99.1M 1s
    ##  70050K .......... .......... .......... .......... .......... 52%  119M 1s
    ##  70100K .......... .......... .......... .......... .......... 52% 89.6M 1s
    ##  70150K .......... .......... .......... .......... .......... 52%  120M 1s
    ##  70200K .......... .......... .......... .......... .......... 52%  101M 1s
    ##  70250K .......... .......... .......... .......... .......... 52%  115M 1s
    ##  70300K .......... .......... .......... .......... .......... 52%  120M 1s
    ##  70350K .......... .......... .......... .......... .......... 52% 94.5M 1s
    ##  70400K .......... .......... .......... .......... .......... 52% 87.9M 1s
    ##  70450K .......... .......... .......... .......... .......... 52% 85.9M 1s
    ##  70500K .......... .......... .......... .......... .......... 52% 86.8M 1s
    ##  70550K .......... .......... .......... .......... .......... 52% 95.1M 1s
    ##  70600K .......... .......... .......... .......... .......... 52%  105M 1s
    ##  70650K .......... .......... .......... .......... .......... 52% 93.2M 1s
    ##  70700K .......... .......... .......... .......... .......... 52%  101M 1s
    ##  70750K .......... .......... .......... .......... .......... 52%  106M 1s
    ##  70800K .......... .......... .......... .......... .......... 52% 84.2M 1s
    ##  70850K .......... .......... .......... .......... .......... 52%  103M 1s
    ##  70900K .......... .......... .......... .......... .......... 52%  108M 1s
    ##  70950K .......... .......... .......... .......... .......... 52% 72.0M 1s
    ##  71000K .......... .......... .......... .......... .......... 52% 84.7M 1s
    ##  71050K .......... .......... .......... .......... .......... 52%  119M 1s
    ##  71100K .......... .......... .......... .......... .......... 52% 68.7M 1s
    ##  71150K .......... .......... .......... .......... .......... 52% 83.3M 1s
    ##  71200K .......... .......... .......... .......... .......... 52% 94.2M 1s
    ##  71250K .......... .......... .......... .......... .......... 52%  101M 1s
    ##  71300K .......... .......... .......... .......... .......... 52%  107M 1s
    ##  71350K .......... .......... .......... .......... .......... 52%  122M 1s
    ##  71400K .......... .......... .......... .......... .......... 53% 92.4M 1s
    ##  71450K .......... .......... .......... .......... .......... 53%  102M 1s
    ##  71500K .......... .......... .......... .......... .......... 53%  102M 1s
    ##  71550K .......... .......... .......... .......... .......... 53% 71.8M 1s
    ##  71600K .......... .......... .......... .......... .......... 53% 99.8M 1s
    ##  71650K .......... .......... .......... .......... .......... 53%  120M 1s
    ##  71700K .......... .......... .......... .......... .......... 53% 49.6M 1s
    ##  71750K .......... .......... .......... .......... .......... 53% 95.9M 1s
    ##  71800K .......... .......... .......... .......... .......... 53% 99.8M 1s
    ##  71850K .......... .......... .......... .......... .......... 53% 82.4M 1s
    ##  71900K .......... .......... .......... .......... .......... 53% 92.2M 1s
    ##  71950K .......... .......... .......... .......... .......... 53%  103M 1s
    ##  72000K .......... .......... .......... .......... .......... 53% 90.4M 1s
    ##  72050K .......... .......... .......... .......... .......... 53% 95.6M 1s
    ##  72100K .......... .......... .......... .......... .......... 53% 99.2M 1s
    ##  72150K .......... .......... .......... .......... .......... 53% 99.9M 1s
    ##  72200K .......... .......... .......... .......... .......... 53% 92.1M 1s
    ##  72250K .......... .......... .......... .......... .......... 53%  111M 1s
    ##  72300K .......... .......... .......... .......... .......... 53% 78.7M 1s
    ##  72350K .......... .......... .......... .......... .......... 53% 78.8M 1s
    ##  72400K .......... .......... .......... .......... .......... 53%  120M 1s
    ##  72450K .......... .......... .......... .......... .......... 53% 89.9M 1s
    ##  72500K .......... .......... .......... .......... .......... 53% 83.6M 1s
    ##  72550K .......... .......... .......... .......... .......... 53%  107M 1s
    ##  72600K .......... .......... .......... .......... .......... 53% 93.5M 1s
    ##  72650K .......... .......... .......... .......... .......... 53% 73.0M 1s
    ##  72700K .......... .......... .......... .......... .......... 53%  108M 1s
    ##  72750K .......... .......... .......... .......... .......... 54%  113M 1s
    ##  72800K .......... .......... .......... .......... .......... 54% 65.1M 1s
    ##  72850K .......... .......... .......... .......... .......... 54%  131M 1s
    ##  72900K .......... .......... .......... .......... .......... 54% 99.5M 1s
    ##  72950K .......... .......... .......... .......... .......... 54% 84.4M 1s
    ##  73000K .......... .......... .......... .......... .......... 54%  121M 1s
    ##  73050K .......... .......... .......... .......... .......... 54%  106M 1s
    ##  73100K .......... .......... .......... .......... .......... 54% 85.8M 1s
    ##  73150K .......... .......... .......... .......... .......... 54%  141M 1s
    ##  73200K .......... .......... .......... .......... .......... 54%  103M 1s
    ##  73250K .......... .......... .......... .......... .......... 54%  103M 1s
    ##  73300K .......... .......... .......... .......... .......... 54%  118M 1s
    ##  73350K .......... .......... .......... .......... .......... 54% 95.7M 1s
    ##  73400K .......... .......... .......... .......... .......... 54%  121M 1s
    ##  73450K .......... .......... .......... .......... .......... 54%  115M 1s
    ##  73500K .......... .......... .......... .......... .......... 54%  120M 1s
    ##  73550K .......... .......... .......... .......... .......... 54%  103M 1s
    ##  73600K .......... .......... .......... .......... .......... 54% 89.3M 1s
    ##  73650K .......... .......... .......... .......... .......... 54%  124M 1s
    ##  73700K .......... .......... .......... .......... .......... 54%  118M 1s
    ##  73750K .......... .......... .......... .......... .......... 54% 93.6M 1s
    ##  73800K .......... .......... .......... .......... .......... 54%  102M 1s
    ##  73850K .......... .......... .......... .......... .......... 54%  138M 1s
    ##  73900K .......... .......... .......... .......... .......... 54% 99.3M 1s
    ##  73950K .......... .......... .......... .......... .......... 54%  105M 1s
    ##  74000K .......... .......... .......... .......... .......... 54%  114M 1s
    ##  74050K .......... .......... .......... .......... .......... 54% 96.8M 1s
    ##  74100K .......... .......... .......... .......... .......... 55%  116M 1s
    ##  74150K .......... .......... .......... .......... .......... 55% 8.53M 1s
    ##  74200K .......... .......... .......... .......... .......... 55% 99.6M 1s
    ##  74250K .......... .......... .......... .......... .......... 55% 30.7M 1s
    ##  74300K .......... .......... .......... .......... .......... 55%  105M 1s
    ##  74350K .......... .......... .......... .......... .......... 55% 85.4M 1s
    ##  74400K .......... .......... .......... .......... .......... 55% 90.0M 1s
    ##  74450K .......... .......... .......... .......... .......... 55%  103M 1s
    ##  74500K .......... .......... .......... .......... .......... 55%  135M 1s
    ##  74550K .......... .......... .......... .......... .......... 55%  139M 1s
    ##  74600K .......... .......... .......... .......... .......... 55% 38.7M 1s
    ##  74650K .......... .......... .......... .......... .......... 55% 35.5M 1s
    ##  74700K .......... .......... .......... .......... .......... 55% 54.2M 1s
    ##  74750K .......... .......... .......... .......... .......... 55% 78.3M 1s
    ##  74800K .......... .......... .......... .......... .......... 55% 38.3M 1s
    ##  74850K .......... .......... .......... .......... .......... 55% 74.1M 1s
    ##  74900K .......... .......... .......... .......... .......... 55% 88.1M 1s
    ##  74950K .......... .......... .......... .......... .......... 55% 83.3M 1s
    ##  75000K .......... .......... .......... .......... .......... 55% 90.0M 1s
    ##  75050K .......... .......... .......... .......... .......... 55%  101M 1s
    ##  75100K .......... .......... .......... .......... .......... 55% 96.3M 1s
    ##  75150K .......... .......... .......... .......... .......... 55% 97.3M 1s
    ##  75200K .......... .......... .......... .......... .......... 55% 73.3M 1s
    ##  75250K .......... .......... .......... .......... .......... 55% 94.9M 1s
    ##  75300K .......... .......... .......... .......... .......... 55% 98.7M 1s
    ##  75350K .......... .......... .......... .......... .......... 55% 8.71M 1s
    ##  75400K .......... .......... .......... .......... .......... 55% 78.4M 1s
    ##  75450K .......... .......... .......... .......... .......... 56% 32.4M 1s
    ##  75500K .......... .......... .......... .......... .......... 56% 78.3M 1s
    ##  75550K .......... .......... .......... .......... .......... 56% 34.4M 1s
    ##  75600K .......... .......... .......... .......... .......... 56%  104M 1s
    ##  75650K .......... .......... .......... .......... .......... 56%  100M 1s
    ##  75700K .......... .......... .......... .......... .......... 56% 79.7M 1s
    ##  75750K .......... .......... .......... .......... .......... 56%  110M 1s
    ##  75800K .......... .......... .......... .......... .......... 56% 29.8M 1s
    ##  75850K .......... .......... .......... .......... .......... 56% 75.6M 1s
    ##  75900K .......... .......... .......... .......... .......... 56% 91.7M 1s
    ##  75950K .......... .......... .......... .......... .......... 56% 94.9M 1s
    ##  76000K .......... .......... .......... .......... .......... 56% 65.4M 1s
    ##  76050K .......... .......... .......... .......... .......... 56% 96.2M 1s
    ##  76100K .......... .......... .......... .......... .......... 56% 72.2M 1s
    ##  76150K .......... .......... .......... .......... .......... 56% 67.2M 1s
    ##  76200K .......... .......... .......... .......... .......... 56% 34.6M 1s
    ##  76250K .......... .......... .......... .......... .......... 56% 51.4M 1s
    ##  76300K .......... .......... .......... .......... .......... 56% 86.5M 1s
    ##  76350K .......... .......... .......... .......... .......... 56% 88.7M 1s
    ##  76400K .......... .......... .......... .......... .......... 56% 92.2M 1s
    ##  76450K .......... .......... .......... .......... .......... 56% 85.4M 1s
    ##  76500K .......... .......... .......... .......... .......... 56% 99.6M 1s
    ##  76550K .......... .......... .......... .......... .......... 56%  102M 1s
    ##  76600K .......... .......... .......... .......... .......... 56% 75.5M 1s
    ##  76650K .......... .......... .......... .......... .......... 56% 86.1M 1s
    ##  76700K .......... .......... .......... .......... .......... 56% 69.3M 1s
    ##  76750K .......... .......... .......... .......... .......... 56% 72.3M 1s
    ##  76800K .......... .......... .......... .......... .......... 57% 88.5M 1s
    ##  76850K .......... .......... .......... .......... .......... 57% 87.3M 1s
    ##  76900K .......... .......... .......... .......... .......... 57% 75.0M 1s
    ##  76950K .......... .......... .......... .......... .......... 57% 98.7M 1s
    ##  77000K .......... .......... .......... .......... .......... 57% 87.0M 1s
    ##  77050K .......... .......... .......... .......... .......... 57% 81.2M 1s
    ##  77100K .......... .......... .......... .......... .......... 57% 83.4M 1s
    ##  77150K .......... .......... .......... .......... .......... 57% 44.0M 1s
    ##  77200K .......... .......... .......... .......... .......... 57% 71.8M 1s
    ##  77250K .......... .......... .......... .......... .......... 57%  104M 1s
    ##  77300K .......... .......... .......... .......... .......... 57% 96.9M 1s
    ##  77350K .......... .......... .......... .......... .......... 57% 44.7M 1s
    ##  77400K .......... .......... .......... .......... .......... 57% 96.1M 1s
    ##  77450K .......... .......... .......... .......... .......... 57%  104M 1s
    ##  77500K .......... .......... .......... .......... .......... 57% 87.1M 1s
    ##  77550K .......... .......... .......... .......... .......... 57%  112M 1s
    ##  77600K .......... .......... .......... .......... .......... 57%  102M 1s
    ##  77650K .......... .......... .......... .......... .......... 57% 63.7M 1s
    ##  77700K .......... .......... .......... .......... .......... 57% 98.3M 1s
    ##  77750K .......... .......... .......... .......... .......... 57% 98.5M 1s
    ##  77800K .......... .......... .......... .......... .......... 57% 63.4M 1s
    ##  77850K .......... .......... .......... .......... .......... 57%  121M 1s
    ##  77900K .......... .......... .......... .......... .......... 57%  106M 1s
    ##  77950K .......... .......... .......... .......... .......... 57%  102M 1s
    ##  78000K .......... .......... .......... .......... .......... 57% 17.1M 1s
    ##  78050K .......... .......... .......... .......... .......... 57% 93.8M 1s
    ##  78100K .......... .......... .......... .......... .......... 58% 89.4M 1s
    ##  78150K .......... .......... .......... .......... .......... 58%  111M 1s
    ##  78200K .......... .......... .......... .......... .......... 58% 79.0M 1s
    ##  78250K .......... .......... .......... .......... .......... 58% 77.3M 1s
    ##  78300K .......... .......... .......... .......... .......... 58% 82.9M 1s
    ##  78350K .......... .......... .......... .......... .......... 58% 84.4M 1s
    ##  78400K .......... .......... .......... .......... .......... 58% 92.7M 1s
    ##  78450K .......... .......... .......... .......... .......... 58% 81.5M 1s
    ##  78500K .......... .......... .......... .......... .......... 58% 98.2M 1s
    ##  78550K .......... .......... .......... .......... .......... 58% 81.3M 1s
    ##  78600K .......... .......... .......... .......... .......... 58% 86.4M 1s
    ##  78650K .......... .......... .......... .......... .......... 58%  101M 1s
    ##  78700K .......... .......... .......... .......... .......... 58% 80.6M 1s
    ##  78750K .......... .......... .......... .......... .......... 58%  102M 1s
    ##  78800K .......... .......... .......... .......... .......... 58% 89.1M 1s
    ##  78850K .......... .......... .......... .......... .......... 58% 87.2M 1s
    ##  78900K .......... .......... .......... .......... .......... 58% 67.9M 1s
    ##  78950K .......... .......... .......... .......... .......... 58% 85.4M 1s
    ##  79000K .......... .......... .......... .......... .......... 58% 89.3M 1s
    ##  79050K .......... .......... .......... .......... .......... 58% 70.8M 1s
    ##  79100K .......... .......... .......... .......... .......... 58% 98.5M 1s
    ##  79150K .......... .......... .......... .......... .......... 58% 88.1M 1s
    ##  79200K .......... .......... .......... .......... .......... 58% 99.8M 1s
    ##  79250K .......... .......... .......... .......... .......... 58% 90.5M 1s
    ##  79300K .......... .......... .......... .......... .......... 58% 81.2M 1s
    ##  79350K .......... .......... .......... .......... .......... 58% 90.4M 1s
    ##  79400K .......... .......... .......... .......... .......... 58% 87.5M 1s
    ##  79450K .......... .......... .......... .......... .......... 59% 87.2M 1s
    ##  79500K .......... .......... .......... .......... .......... 59% 93.8M 1s
    ##  79550K .......... .......... .......... .......... .......... 59% 96.7M 1s
    ##  79600K .......... .......... .......... .......... .......... 59% 77.6M 1s
    ##  79650K .......... .......... .......... .......... .......... 59% 91.3M 1s
    ##  79700K .......... .......... .......... .......... .......... 59%  102M 1s
    ##  79750K .......... .......... .......... .......... .......... 59% 84.4M 1s
    ##  79800K .......... .......... .......... .......... .......... 59%  111M 1s
    ##  79850K .......... .......... .......... .......... .......... 59% 92.2M 1s
    ##  79900K .......... .......... .......... .......... .......... 59% 91.2M 1s
    ##  79950K .......... .......... .......... .......... .......... 59% 83.9M 1s
    ##  80000K .......... .......... .......... .......... .......... 59% 98.7M 1s
    ##  80050K .......... .......... .......... .......... .......... 59% 90.0M 1s
    ##  80100K .......... .......... .......... .......... .......... 59%  105M 1s
    ##  80150K .......... .......... .......... .......... .......... 59% 91.0M 1s
    ##  80200K .......... .......... .......... .......... .......... 59% 88.1M 1s
    ##  80250K .......... .......... .......... .......... .......... 59%  107M 1s
    ##  80300K .......... .......... .......... .......... .......... 59% 88.6M 1s
    ##  80350K .......... .......... .......... .......... .......... 59%  114M 1s
    ##  80400K .......... .......... .......... .......... .......... 59%  106M 1s
    ##  80450K .......... .......... .......... .......... .......... 59% 88.4M 1s
    ##  80500K .......... .......... .......... .......... .......... 59% 91.0M 1s
    ##  80550K .......... .......... .......... .......... .......... 59% 99.3M 1s
    ##  80600K .......... .......... .......... .......... .......... 59%  111M 1s
    ##  80650K .......... .......... .......... .......... .......... 59% 92.0M 1s
    ##  80700K .......... .......... .......... .......... .......... 59% 93.2M 1s
    ##  80750K .......... .......... .......... .......... .......... 59% 88.7M 1s
    ##  80800K .......... .......... .......... .......... .......... 60%  102M 1s
    ##  80850K .......... .......... .......... .......... .......... 60%  125M 1s
    ##  80900K .......... .......... .......... .......... .......... 60%  109M 1s
    ##  80950K .......... .......... .......... .......... .......... 60%  100M 1s
    ##  81000K .......... .......... .......... .......... .......... 60%  114M 1s
    ##  81050K .......... .......... .......... .......... .......... 60%  104M 1s
    ##  81100K .......... .......... .......... .......... .......... 60%  109M 1s
    ##  81150K .......... .......... .......... .......... .......... 60% 90.3M 1s
    ##  81200K .......... .......... .......... .......... .......... 60% 94.5M 1s
    ##  81250K .......... .......... .......... .......... .......... 60%  105M 1s
    ##  81300K .......... .......... .......... .......... .......... 60% 98.4M 1s
    ##  81350K .......... .......... .......... .......... .......... 60%  124M 1s
    ##  81400K .......... .......... .......... .......... .......... 60% 88.3M 1s
    ##  81450K .......... .......... .......... .......... .......... 60%  101M 1s
    ##  81500K .......... .......... .......... .......... .......... 60% 88.5M 1s
    ##  81550K .......... .......... .......... .......... .......... 60%  106M 1s
    ##  81600K .......... .......... .......... .......... .......... 60%  100M 1s
    ##  81650K .......... .......... .......... .......... .......... 60% 81.1M 1s
    ##  81700K .......... .......... .......... .......... .......... 60%  105M 1s
    ##  81750K .......... .......... .......... .......... .......... 60% 96.5M 1s
    ##  81800K .......... .......... .......... .......... .......... 60%  113M 1s
    ##  81850K .......... .......... .......... .......... .......... 60%  120M 1s
    ##  81900K .......... .......... .......... .......... .......... 60%  102M 1s
    ##  81950K .......... .......... .......... .......... .......... 60% 98.9M 1s
    ##  82000K .......... .......... .......... .......... .......... 60%  102M 1s
    ##  82050K .......... .......... .......... .......... .......... 60%  118M 1s
    ##  82100K .......... .......... .......... .......... .......... 60% 96.7M 1s
    ##  82150K .......... .......... .......... .......... .......... 61% 96.5M 1s
    ##  82200K .......... .......... .......... .......... .......... 61%  101M 1s
    ##  82250K .......... .......... .......... .......... .......... 61%  104M 1s
    ##  82300K .......... .......... .......... .......... .......... 61%  110M 1s
    ##  82350K .......... .......... .......... .......... .......... 61%  124M 1s
    ##  82400K .......... .......... .......... .......... .......... 61% 99.3M 1s
    ##  82450K .......... .......... .......... .......... .......... 61%  124M 1s
    ##  82500K .......... .......... .......... .......... .......... 61%  106M 1s
    ##  82550K .......... .......... .......... .......... .......... 61%  109M 1s
    ##  82600K .......... .......... .......... .......... .......... 61% 85.9M 1s
    ##  82650K .......... .......... .......... .......... .......... 61%  123M 1s
    ##  82700K .......... .......... .......... .......... .......... 61% 98.6M 1s
    ##  82750K .......... .......... .......... .......... .......... 61%  106M 1s
    ##  82800K .......... .......... .......... .......... .......... 61%  117M 1s
    ##  82850K .......... .......... .......... .......... .......... 61%  104M 1s
    ##  82900K .......... .......... .......... .......... .......... 61%  104M 1s
    ##  82950K .......... .......... .......... .......... .......... 61%  113M 1s
    ##  83000K .......... .......... .......... .......... .......... 61%  115M 1s
    ##  83050K .......... .......... .......... .......... .......... 61%  126M 1s
    ##  83100K .......... .......... .......... .......... .......... 61%  117M 1s
    ##  83150K .......... .......... .......... .......... .......... 61%  108M 1s
    ##  83200K .......... .......... .......... .......... .......... 61%  104M 1s
    ##  83250K .......... .......... .......... .......... .......... 61%  118M 1s
    ##  83300K .......... .......... .......... .......... .......... 61%  132M 1s
    ##  83350K .......... .......... .......... .......... .......... 61%  116M 1s
    ##  83400K .......... .......... .......... .......... .......... 61% 55.3M 1s
    ##  83450K .......... .......... .......... .......... .......... 61%  112M 1s
    ##  83500K .......... .......... .......... .......... .......... 62%  136M 1s
    ##  83550K .......... .......... .......... .......... .......... 62%  144M 1s
    ##  83600K .......... .......... .......... .......... .......... 62% 63.9M 1s
    ##  83650K .......... .......... .......... .......... .......... 62%  142M 1s
    ##  83700K .......... .......... .......... .......... .......... 62%  137M 1s
    ##  83750K .......... .......... .......... .......... .......... 62%  116M 1s
    ##  83800K .......... .......... .......... .......... .......... 62%  102M 1s
    ##  83850K .......... .......... .......... .......... .......... 62%  125M 1s
    ##  83900K .......... .......... .......... .......... .......... 62%  126M 1s
    ##  83950K .......... .......... .......... .......... .......... 62% 46.5M 1s
    ##  84000K .......... .......... .......... .......... .......... 62%  101M 1s
    ##  84050K .......... .......... .......... .......... .......... 62%  110M 1s
    ##  84100K .......... .......... .......... .......... .......... 62% 3.75M 1s
    ##  84150K .......... .......... .......... .......... .......... 62%  121M 1s
    ##  84200K .......... .......... .......... .......... .......... 62%  153M 1s
    ##  84250K .......... .......... .......... .......... .......... 62%  142M 1s
    ##  84300K .......... .......... .......... .......... .......... 62%  154M 1s
    ##  84350K .......... .......... .......... .......... .......... 62%  165M 1s
    ##  84400K .......... .......... .......... .......... .......... 62%  143M 1s
    ##  84450K .......... .......... .......... .......... .......... 62%  144M 1s
    ##  84500K .......... .......... .......... .......... .......... 62%  150M 1s
    ##  84550K .......... .......... .......... .......... .......... 62%  146M 1s
    ##  84600K .......... .......... .......... .......... .......... 62%  151M 1s
    ##  84650K .......... .......... .......... .......... .......... 62%  133M 1s
    ##  84700K .......... .......... .......... .......... .......... 62%  133M 1s
    ##  84750K .......... .......... .......... .......... .......... 62%  132M 1s
    ##  84800K .......... .......... .......... .......... .......... 62%  152M 1s
    ##  84850K .......... .......... .......... .......... .......... 63%  142M 1s
    ##  84900K .......... .......... .......... .......... .......... 63%  152M 1s
    ##  84950K .......... .......... .......... .......... .......... 63%  163M 1s
    ##  85000K .......... .......... .......... .......... .......... 63%  106M 1s
    ##  85050K .......... .......... .......... .......... .......... 63%  101M 1s
    ##  85100K .......... .......... .......... .......... .......... 63%  156M 1s
    ##  85150K .......... .......... .......... .......... .......... 63%  146M 1s
    ##  85200K .......... .......... .......... .......... .......... 63%  146M 1s
    ##  85250K .......... .......... .......... .......... .......... 63%  157M 1s
    ##  85300K .......... .......... .......... .......... .......... 63% 13.0M 1s
    ##  85350K .......... .......... .......... .......... .......... 63%  155M 1s
    ##  85400K .......... .......... .......... .......... .......... 63%  170M 1s
    ##  85450K .......... .......... .......... .......... .......... 63%  145M 1s
    ##  85500K .......... .......... .......... .......... .......... 63%  156M 1s
    ##  85550K .......... .......... .......... .......... .......... 63% 31.4M 1s
    ##  85600K .......... .......... .......... .......... .......... 63% 51.4M 1s
    ##  85650K .......... .......... .......... .......... .......... 63% 42.7M 1s
    ##  85700K .......... .......... .......... .......... .......... 63%  120M 1s
    ##  85750K .......... .......... .......... .......... .......... 63% 72.4M 1s
    ##  85800K .......... .......... .......... .......... .......... 63%  140M 1s
    ##  85850K .......... .......... .......... .......... .......... 63% 80.3M 1s
    ##  85900K .......... .......... .......... .......... .......... 63% 57.4M 1s
    ##  85950K .......... .......... .......... .......... .......... 63% 56.3M 1s
    ##  86000K .......... .......... .......... .......... .......... 63%  110M 1s
    ##  86050K .......... .......... .......... .......... .......... 63% 34.0M 1s
    ##  86100K .......... .......... .......... .......... .......... 63% 96.7M 1s
    ##  86150K .......... .......... .......... .......... .......... 63% 50.0M 1s
    ##  86200K .......... .......... .......... .......... .......... 64%  127M 1s
    ##  86250K .......... .......... .......... .......... .......... 64%  160M 1s
    ##  86300K .......... .......... .......... .......... .......... 64%  150M 1s
    ##  86350K .......... .......... .......... .......... .......... 64% 43.6M 1s
    ##  86400K .......... .......... .......... .......... .......... 64%  120M 1s
    ##  86450K .......... .......... .......... .......... .......... 64%  129M 1s
    ##  86500K .......... .......... .......... .......... .......... 64% 94.5M 1s
    ##  86550K .......... .......... .......... .......... .......... 64% 42.0M 1s
    ##  86600K .......... .......... .......... .......... .......... 64% 26.4M 1s
    ##  86650K .......... .......... .......... .......... .......... 64%  121M 1s
    ##  86700K .......... .......... .......... .......... .......... 64%  145M 1s
    ##  86750K .......... .......... .......... .......... .......... 64%  165M 1s
    ##  86800K .......... .......... .......... .......... .......... 64%  159M 1s
    ##  86850K .......... .......... .......... .......... .......... 64%  156M 1s
    ##  86900K .......... .......... .......... .......... .......... 64% 46.4M 1s
    ##  86950K .......... .......... .......... .......... .......... 64% 19.3M 1s
    ##  87000K .......... .......... .......... .......... .......... 64%  137M 1s
    ##  87050K .......... .......... .......... .......... .......... 64% 25.3M 1s
    ##  87100K .......... .......... .......... .......... .......... 64%  155M 1s
    ##  87150K .......... .......... .......... .......... .......... 64%  160M 1s
    ##  87200K .......... .......... .......... .......... .......... 64%  130M 1s
    ##  87250K .......... .......... .......... .......... .......... 64% 21.8M 1s
    ##  87300K .......... .......... .......... .......... .......... 64%  129M 1s
    ##  87350K .......... .......... .......... .......... .......... 64% 39.6M 1s
    ##  87400K .......... .......... .......... .......... .......... 64% 55.1M 1s
    ##  87450K .......... .......... .......... .......... .......... 64%  128M 1s
    ##  87500K .......... .......... .......... .......... .......... 64%  130M 1s
    ##  87550K .......... .......... .......... .......... .......... 65% 52.6M 1s
    ##  87600K .......... .......... .......... .......... .......... 65%  109M 1s
    ##  87650K .......... .......... .......... .......... .......... 65% 99.4M 1s
    ##  87700K .......... .......... .......... .......... .......... 65% 47.9M 1s
    ##  87750K .......... .......... .......... .......... .......... 65% 37.4M 1s
    ##  87800K .......... .......... .......... .......... .......... 65% 61.0M 1s
    ##  87850K .......... .......... .......... .......... .......... 65% 31.2M 1s
    ##  87900K .......... .......... .......... .......... .......... 65% 85.6M 1s
    ##  87950K .......... .......... .......... .......... .......... 65%  112M 1s
    ##  88000K .......... .......... .......... .......... .......... 65% 56.1M 1s
    ##  88050K .......... .......... .......... .......... .......... 65% 54.0M 1s
    ##  88100K .......... .......... .......... .......... .......... 65% 47.0M 1s
    ##  88150K .......... .......... .......... .......... .......... 65% 82.9M 1s
    ##  88200K .......... .......... .......... .......... .......... 65% 37.3M 1s
    ##  88250K .......... .......... .......... .......... .......... 65%  114M 1s
    ##  88300K .......... .......... .......... .......... .......... 65% 28.8M 1s
    ##  88350K .......... .......... .......... .......... .......... 65%  119M 1s
    ##  88400K .......... .......... .......... .......... .......... 65%  152M 1s
    ##  88450K .......... .......... .......... .......... .......... 65% 64.4M 1s
    ##  88500K .......... .......... .......... .......... .......... 65% 48.8M 1s
    ##  88550K .......... .......... .......... .......... .......... 65% 56.9M 1s
    ##  88600K .......... .......... .......... .......... .......... 65%  103M 1s
    ##  88650K .......... .......... .......... .......... .......... 65% 36.7M 1s
    ##  88700K .......... .......... .......... .......... .......... 65%  117M 1s
    ##  88750K .......... .......... .......... .......... .......... 65% 66.9M 1s
    ##  88800K .......... .......... .......... .......... .......... 65%  142M 1s
    ##  88850K .......... .......... .......... .......... .......... 65% 63.1M 1s
    ##  88900K .......... .......... .......... .......... .......... 66% 97.9M 1s
    ##  88950K .......... .......... .......... .......... .......... 66% 80.1M 1s
    ##  89000K .......... .......... .......... .......... .......... 66% 59.1M 1s
    ##  89050K .......... .......... .......... .......... .......... 66%  116M 1s
    ##  89100K .......... .......... .......... .......... .......... 66% 82.5M 1s
    ##  89150K .......... .......... .......... .......... .......... 66% 89.2M 1s
    ##  89200K .......... .......... .......... .......... .......... 66% 87.9M 1s
    ##  89250K .......... .......... .......... .......... .......... 66% 92.1M 1s
    ##  89300K .......... .......... .......... .......... .......... 66%  101M 1s
    ##  89350K .......... .......... .......... .......... .......... 66% 66.9M 1s
    ##  89400K .......... .......... .......... .......... .......... 66% 81.3M 1s
    ##  89450K .......... .......... .......... .......... .......... 66%  162M 1s
    ##  89500K .......... .......... .......... .......... .......... 66% 85.7M 1s
    ##  89550K .......... .......... .......... .......... .......... 66% 40.5M 1s
    ##  89600K .......... .......... .......... .......... .......... 66%  118M 1s
    ##  89650K .......... .......... .......... .......... .......... 66% 77.3M 1s
    ##  89700K .......... .......... .......... .......... .......... 66% 63.3M 1s
    ##  89750K .......... .......... .......... .......... .......... 66%  147M 1s
    ##  89800K .......... .......... .......... .......... .......... 66% 52.9M 1s
    ##  89850K .......... .......... .......... .......... .......... 66% 34.4M 1s
    ##  89900K .......... .......... .......... .......... .......... 66% 26.9M 1s
    ##  89950K .......... .......... .......... .......... .......... 66%  117M 1s
    ##  90000K .......... .......... .......... .......... .......... 66%  141M 1s
    ##  90050K .......... .......... .......... .......... .......... 66%  187M 1s
    ##  90100K .......... .......... .......... .......... .......... 66%  135M 1s
    ##  90150K .......... .......... .......... .......... .......... 66%  161M 1s
    ##  90200K .......... .......... .......... .......... .......... 66% 24.1M 1s
    ##  90250K .......... .......... .......... .......... .......... 67% 88.1M 1s
    ##  90300K .......... .......... .......... .......... .......... 67% 90.4M 1s
    ##  90350K .......... .......... .......... .......... .......... 67% 46.2M 1s
    ##  90400K .......... .......... .......... .......... .......... 67%  109M 1s
    ##  90450K .......... .......... .......... .......... .......... 67%  163M 1s
    ##  90500K .......... .......... .......... .......... .......... 67%  117M 1s
    ##  90550K .......... .......... .......... .......... .......... 67% 57.6M 1s
    ##  90600K .......... .......... .......... .......... .......... 67%  128M 1s
    ##  90650K .......... .......... .......... .......... .......... 67% 34.5M 1s
    ##  90700K .......... .......... .......... .......... .......... 67%  126M 1s
    ##  90750K .......... .......... .......... .......... .......... 67%  147M 1s
    ##  90800K .......... .......... .......... .......... .......... 67%  142M 1s
    ##  90850K .......... .......... .......... .......... .......... 67% 30.8M 1s
    ##  90900K .......... .......... .......... .......... .......... 67%  117M 1s
    ##  90950K .......... .......... .......... .......... .......... 67%  157M 1s
    ##  91000K .......... .......... .......... .......... .......... 67%  142M 1s
    ##  91050K .......... .......... .......... .......... .......... 67%  174M 1s
    ##  91100K .......... .......... .......... .......... .......... 67%  130M 1s
    ##  91150K .......... .......... .......... .......... .......... 67% 47.7M 1s
    ##  91200K .......... .......... .......... .......... .......... 67% 41.4M 1s
    ##  91250K .......... .......... .......... .......... .......... 67% 74.0M 1s
    ##  91300K .......... .......... .......... .......... .......... 67% 42.5M 1s
    ##  91350K .......... .......... .......... .......... .......... 67% 95.2M 1s
    ##  91400K .......... .......... .......... .......... .......... 67% 56.2M 1s
    ##  91450K .......... .......... .......... .......... .......... 67%  126M 1s
    ##  91500K .......... .......... .......... .......... .......... 67% 54.1M 1s
    ##  91550K .......... .......... .......... .......... .......... 67%  127M 1s
    ##  91600K .......... .......... .......... .......... .......... 68% 52.3M 1s
    ##  91650K .......... .......... .......... .......... .......... 68%  133M 1s
    ##  91700K .......... .......... .......... .......... .......... 68% 44.2M 1s
    ##  91750K .......... .......... .......... .......... .......... 68% 95.4M 1s
    ##  91800K .......... .......... .......... .......... .......... 68% 73.6M 1s
    ##  91850K .......... .......... .......... .......... .......... 68% 47.4M 1s
    ##  91900K .......... .......... .......... .......... .......... 68% 99.3M 1s
    ##  91950K .......... .......... .......... .......... .......... 68% 68.5M 1s
    ##  92000K .......... .......... .......... .......... .......... 68% 50.1M 1s
    ##  92050K .......... .......... .......... .......... .......... 68%  131M 1s
    ##  92100K .......... .......... .......... .......... .......... 68% 97.8M 1s
    ##  92150K .......... .......... .......... .......... .......... 68% 73.5M 1s
    ##  92200K .......... .......... .......... .......... .......... 68% 62.9M 1s
    ##  92250K .......... .......... .......... .......... .......... 68% 68.9M 1s
    ##  92300K .......... .......... .......... .......... .......... 68% 54.1M 1s
    ##  92350K .......... .......... .......... .......... .......... 68% 79.5M 1s
    ##  92400K .......... .......... .......... .......... .......... 68%  102M 1s
    ##  92450K .......... .......... .......... .......... .......... 68% 73.2M 1s
    ##  92500K .......... .......... .......... .......... .......... 68% 62.7M 1s
    ##  92550K .......... .......... .......... .......... .......... 68% 70.8M 1s
    ##  92600K .......... .......... .......... .......... .......... 68%  101M 1s
    ##  92650K .......... .......... .......... .......... .......... 68% 66.5M 1s
    ##  92700K .......... .......... .......... .......... .......... 68% 64.8M 1s
    ##  92750K .......... .......... .......... .......... .......... 68% 95.3M 1s
    ##  92800K .......... .......... .......... .......... .......... 68% 90.8M 1s
    ##  92850K .......... .......... .......... .......... .......... 68% 74.5M 1s
    ##  92900K .......... .......... .......... .......... .......... 68% 67.1M 1s
    ##  92950K .......... .......... .......... .......... .......... 69% 64.7M 1s
    ##  93000K .......... .......... .......... .......... .......... 69% 81.2M 1s
    ##  93050K .......... .......... .......... .......... .......... 69%  142M 1s
    ##  93100K .......... .......... .......... .......... .......... 69% 81.0M 1s
    ##  93150K .......... .......... .......... .......... .......... 69% 65.5M 1s
    ##  93200K .......... .......... .......... .......... .......... 69% 68.5M 1s
    ##  93250K .......... .......... .......... .......... .......... 69% 76.6M 1s
    ##  93300K .......... .......... .......... .......... .......... 69%  119M 1s
    ##  93350K .......... .......... .......... .......... .......... 69% 69.3M 1s
    ##  93400K .......... .......... .......... .......... .......... 69% 73.4M 1s
    ##  93450K .......... .......... .......... .......... .......... 69% 60.5M 1s
    ##  93500K .......... .......... .......... .......... .......... 69% 76.9M 1s
    ##  93550K .......... .......... .......... .......... .......... 69%  124M 1s
    ##  93600K .......... .......... .......... .......... .......... 69% 69.6M 1s
    ##  93650K .......... .......... .......... .......... .......... 69% 88.6M 1s
    ##  93700K .......... .......... .......... .......... .......... 69% 81.0M 1s
    ##  93750K .......... .......... .......... .......... .......... 69% 67.8M 1s
    ##  93800K .......... .......... .......... .......... .......... 69%  119M 1s
    ##  93850K .......... .......... .......... .......... .......... 69% 86.6M 1s
    ##  93900K .......... .......... .......... .......... .......... 69% 58.2M 1s
    ##  93950K .......... .......... .......... .......... .......... 69%  115M 1s
    ##  94000K .......... .......... .......... .......... .......... 69%  110M 1s
    ##  94050K .......... .......... .......... .......... .......... 69% 86.0M 1s
    ##  94100K .......... .......... .......... .......... .......... 69% 65.6M 1s
    ##  94150K .......... .......... .......... .......... .......... 69%  111M 1s
    ##  94200K .......... .......... .......... .......... .......... 69% 67.3M 1s
    ##  94250K .......... .......... .......... .......... .......... 69%  149M 1s
    ##  94300K .......... .......... .......... .......... .......... 70% 60.5M 1s
    ##  94350K .......... .......... .......... .......... .......... 70% 67.3M 1s
    ##  94400K .......... .......... .......... .......... .......... 70% 89.2M 1s
    ##  94450K .......... .......... .......... .......... .......... 70%  101M 1s
    ##  94500K .......... .......... .......... .......... .......... 70% 79.7M 1s
    ##  94550K .......... .......... .......... .......... .......... 70% 83.1M 1s
    ##  94600K .......... .......... .......... .......... .......... 70% 92.2M 1s
    ##  94650K .......... .......... .......... .......... .......... 70% 75.3M 1s
    ##  94700K .......... .......... .......... .......... .......... 70%  134M 1s
    ##  94750K .......... .......... .......... .......... .......... 70% 82.8M 1s
    ##  94800K .......... .......... .......... .......... .......... 70% 19.1M 1s
    ##  94850K .......... .......... .......... .......... .......... 70%  110M 1s
    ##  94900K .......... .......... .......... .......... .......... 70%  114M 1s
    ##  94950K .......... .......... .......... .......... .......... 70%  172M 1s
    ##  95000K .......... .......... .......... .......... .......... 70%  159M 1s
    ##  95050K .......... .......... .......... .......... .......... 70%  195M 1s
    ##  95100K .......... .......... .......... .......... .......... 70% 28.8M 1s
    ##  95150K .......... .......... .......... .......... .......... 70% 41.7M 1s
    ##  95200K .......... .......... .......... .......... .......... 70% 43.6M 1s
    ##  95250K .......... .......... .......... .......... .......... 70%  157M 1s
    ##  95300K .......... .......... .......... .......... .......... 70% 49.1M 1s
    ##  95350K .......... .......... .......... .......... .......... 70%  125M 1s
    ##  95400K .......... .......... .......... .......... .......... 70%  140M 1s
    ##  95450K .......... .......... .......... .......... .......... 70%  152M 1s
    ##  95500K .......... .......... .......... .......... .......... 70% 11.7M 1s
    ##  95550K .......... .......... .......... .......... .......... 70%  170M 1s
    ##  95600K .......... .......... .......... .......... .......... 70%  151M 1s
    ##  95650K .......... .......... .......... .......... .......... 71%  182M 1s
    ##  95700K .......... .......... .......... .......... .......... 71%  149M 1s
    ##  95750K .......... .......... .......... .......... .......... 71%  163M 1s
    ##  95800K .......... .......... .......... .......... .......... 71% 36.5M 1s
    ##  95850K .......... .......... .......... .......... .......... 71% 27.7M 1s
    ##  95900K .......... .......... .......... .......... .......... 71% 89.4M 1s
    ##  95950K .......... .......... .......... .......... .......... 71% 41.7M 1s
    ##  96000K .......... .......... .......... .......... .......... 71% 59.7M 1s
    ##  96050K .......... .......... .......... .......... .......... 71% 73.1M 1s
    ##  96100K .......... .......... .......... .......... .......... 71% 86.0M 1s
    ##  96150K .......... .......... .......... .......... .......... 71% 97.7M 1s
    ##  96200K .......... .......... .......... .......... .......... 71% 92.0M 1s
    ##  96250K .......... .......... .......... .......... .......... 71% 87.0M 1s
    ##  96300K .......... .......... .......... .......... .......... 71% 66.9M 1s
    ##  96350K .......... .......... .......... .......... .......... 71% 43.3M 1s
    ##  96400K .......... .......... .......... .......... .......... 71% 75.7M 1s
    ##  96450K .......... .......... .......... .......... .......... 71% 82.7M 1s
    ##  96500K .......... .......... .......... .......... .......... 71% 99.0M 1s
    ##  96550K .......... .......... .......... .......... .......... 71% 50.5M 1s
    ##  96600K .......... .......... .......... .......... .......... 71% 33.9M 1s
    ##  96650K .......... .......... .......... .......... .......... 71% 83.8M 1s
    ##  96700K .......... .......... .......... .......... .......... 71% 90.0M 1s
    ##  96750K .......... .......... .......... .......... .......... 71% 81.5M 1s
    ##  96800K .......... .......... .......... .......... .......... 71% 92.8M 1s
    ##  96850K .......... .......... .......... .......... .......... 71%  114M 1s
    ##  96900K .......... .......... .......... .......... .......... 71% 79.8M 1s
    ##  96950K .......... .......... .......... .......... .......... 71%  103M 1s
    ##  97000K .......... .......... .......... .......... .......... 72% 74.3M 1s
    ##  97050K .......... .......... .......... .......... .......... 72% 95.9M 1s
    ##  97100K .......... .......... .......... .......... .......... 72% 38.1M 1s
    ##  97150K .......... .......... .......... .......... .......... 72% 97.6M 1s
    ##  97200K .......... .......... .......... .......... .......... 72% 82.7M 1s
    ##  97250K .......... .......... .......... .......... .......... 72% 88.5M 1s
    ##  97300K .......... .......... .......... .......... .......... 72%  111M 1s
    ##  97350K .......... .......... .......... .......... .......... 72%  120M 1s
    ##  97400K .......... .......... .......... .......... .......... 72% 46.9M 1s
    ##  97450K .......... .......... .......... .......... .......... 72% 52.1M 1s
    ##  97500K .......... .......... .......... .......... .......... 72% 60.3M 1s
    ##  97550K .......... .......... .......... .......... .......... 72% 81.8M 1s
    ##  97600K .......... .......... .......... .......... .......... 72% 53.4M 1s
    ##  97650K .......... .......... .......... .......... .......... 72%  102M 1s
    ##  97700K .......... .......... .......... .......... .......... 72%  112M 1s
    ##  97750K .......... .......... .......... .......... .......... 72% 24.0M 1s
    ##  97800K .......... .......... .......... .......... .......... 72% 83.9M 1s
    ##  97850K .......... .......... .......... .......... .......... 72% 33.2M 1s
    ##  97900K .......... .......... .......... .......... .......... 72% 87.0M 1s
    ##  97950K .......... .......... .......... .......... .......... 72%  103M 1s
    ##  98000K .......... .......... .......... .......... .......... 72% 74.1M 1s
    ##  98050K .......... .......... .......... .......... .......... 72%  113M 1s
    ##  98100K .......... .......... .......... .......... .......... 72%  113M 1s
    ##  98150K .......... .......... .......... .......... .......... 72% 70.3M 1s
    ##  98200K .......... .......... .......... .......... .......... 72% 74.1M 1s
    ##  98250K .......... .......... .......... .......... .......... 72% 56.1M 1s
    ##  98300K .......... .......... .......... .......... .......... 72% 67.1M 1s
    ##  98350K .......... .......... .......... .......... .......... 73% 83.0M 1s
    ##  98400K .......... .......... .......... .......... .......... 73% 81.1M 1s
    ##  98450K .......... .......... .......... .......... .......... 73%  108M 1s
    ##  98500K .......... .......... .......... .......... .......... 73% 85.1M 1s
    ##  98550K .......... .......... .......... .......... .......... 73% 47.8M 1s
    ##  98600K .......... .......... .......... .......... .......... 73% 96.5M 1s
    ##  98650K .......... .......... .......... .......... .......... 73% 98.9M 1s
    ##  98700K .......... .......... .......... .......... .......... 73% 69.7M 1s
    ##  98750K .......... .......... .......... .......... .......... 73% 55.7M 1s
    ##  98800K .......... .......... .......... .......... .......... 73% 11.8M 1s
    ##  98850K .......... .......... .......... .......... .......... 73%  107M 1s
    ##  98900K .......... .......... .......... .......... .......... 73%  113M 1s
    ##  98950K .......... .......... .......... .......... .......... 73%  108M 1s
    ##  99000K .......... .......... .......... .......... .......... 73% 88.8M 1s
    ##  99050K .......... .......... .......... .......... .......... 73%  122M 1s
    ##  99100K .......... .......... .......... .......... .......... 73%  104M 1s
    ##  99150K .......... .......... .......... .......... .......... 73% 33.8M 1s
    ##  99200K .......... .......... .......... .......... .......... 73%  119M 1s
    ##  99250K .......... .......... .......... .......... .......... 73% 29.1M 1s
    ##  99300K .......... .......... .......... .......... .......... 73% 84.3M 1s
    ##  99350K .......... .......... .......... .......... .......... 73% 88.5M 1s
    ##  99400K .......... .......... .......... .......... .......... 73%  127M 1s
    ##  99450K .......... .......... .......... .......... .......... 73% 72.6M 1s
    ##  99500K .......... .......... .......... .......... .......... 73%  136M 1s
    ##  99550K .......... .......... .......... .......... .......... 73% 85.8M 1s
    ##  99600K .......... .......... .......... .......... .......... 73% 92.0M 1s
    ##  99650K .......... .......... .......... .......... .......... 73% 29.9M 1s
    ##  99700K .......... .......... .......... .......... .......... 74% 70.4M 1s
    ##  99750K .......... .......... .......... .......... .......... 74% 49.4M 1s
    ##  99800K .......... .......... .......... .......... .......... 74% 86.7M 1s
    ##  99850K .......... .......... .......... .......... .......... 74% 91.0M 1s
    ##  99900K .......... .......... .......... .......... .......... 74%  102M 1s
    ##  99950K .......... .......... .......... .......... .......... 74% 80.2M 1s
    ## 100000K .......... .......... .......... .......... .......... 74% 86.7M 1s
    ## 100050K .......... .......... .......... .......... .......... 74% 78.7M 1s
    ## 100100K .......... .......... .......... .......... .......... 74% 94.4M 1s
    ## 100150K .......... .......... .......... .......... .......... 74% 92.1M 1s
    ## 100200K .......... .......... .......... .......... .......... 74% 45.7M 1s
    ## 100250K .......... .......... .......... .......... .......... 74% 89.5M 1s
    ## 100300K .......... .......... .......... .......... .......... 74% 45.2M 1s
    ## 100350K .......... .......... .......... .......... .......... 74% 95.3M 1s
    ## 100400K .......... .......... .......... .......... .......... 74%  101M 1s
    ## 100450K .......... .......... .......... .......... .......... 74% 65.8M 1s
    ## 100500K .......... .......... .......... .......... .......... 74% 92.5M 1s
    ## 100550K .......... .......... .......... .......... .......... 74%  103M 1s
    ## 100600K .......... .......... .......... .......... .......... 74%  110M 1s
    ## 100650K .......... .......... .......... .......... .......... 74% 46.6M 1s
    ## 100700K .......... .......... .......... .......... .......... 74% 64.5M 1s
    ## 100750K .......... .......... .......... .......... .......... 74% 65.5M 1s
    ## 100800K .......... .......... .......... .......... .......... 74% 99.0M 1s
    ## 100850K .......... .......... .......... .......... .......... 74% 92.5M 1s
    ## 100900K .......... .......... .......... .......... .......... 74% 61.4M 1s
    ## 100950K .......... .......... .......... .......... .......... 74% 15.9M 1s
    ## 101000K .......... .......... .......... .......... .......... 74%  106M 1s
    ## 101050K .......... .......... .......... .......... .......... 75%  107M 1s
    ## 101100K .......... .......... .......... .......... .......... 75%  102M 1s
    ## 101150K .......... .......... .......... .......... .......... 75%  126M 1s
    ## 101200K .......... .......... .......... .......... .......... 75%  115M 1s
    ## 101250K .......... .......... .......... .......... .......... 75% 68.2M 1s
    ## 101300K .......... .......... .......... .......... .......... 75% 36.4M 0s
    ## 101350K .......... .......... .......... .......... .......... 75% 69.8M 0s
    ## 101400K .......... .......... .......... .......... .......... 75% 44.9M 0s
    ## 101450K .......... .......... .......... .......... .......... 75% 87.7M 0s
    ## 101500K .......... .......... .......... .......... .......... 75% 82.9M 0s
    ## 101550K .......... .......... .......... .......... .......... 75%  126M 0s
    ## 101600K .......... .......... .......... .......... .......... 75% 77.9M 0s
    ## 101650K .......... .......... .......... .......... .......... 75%  105M 0s
    ## 101700K .......... .......... .......... .......... .......... 75% 87.1M 0s
    ## 101750K .......... .......... .......... .......... .......... 75% 71.4M 0s
    ## 101800K .......... .......... .......... .......... .......... 75% 92.4M 0s
    ## 101850K .......... .......... .......... .......... .......... 75% 69.2M 0s
    ## 101900K .......... .......... .......... .......... .......... 75% 89.7M 0s
    ## 101950K .......... .......... .......... .......... .......... 75% 65.6M 0s
    ## 102000K .......... .......... .......... .......... .......... 75% 95.8M 0s
    ## 102050K .......... .......... .......... .......... .......... 75% 12.4M 0s
    ## 102100K .......... .......... .......... .......... .......... 75%  108M 0s
    ## 102150K .......... .......... .......... .......... .......... 75%  127M 0s
    ## 102200K .......... .......... .......... .......... .......... 75%  115M 0s
    ## 102250K .......... .......... .......... .......... .......... 75% 97.8M 0s
    ## 102300K .......... .......... .......... .......... .......... 75%  127M 0s
    ## 102350K .......... .......... .......... .......... .......... 75%  131M 0s
    ## 102400K .......... .......... .......... .......... .......... 76% 52.1M 0s
    ## 102450K .......... .......... .......... .......... .......... 76% 61.5M 0s
    ## 102500K .......... .......... .......... .......... .......... 76% 42.4M 0s
    ## 102550K .......... .......... .......... .......... .......... 76% 73.7M 0s
    ## 102600K .......... .......... .......... .......... .......... 76% 95.1M 0s
    ## 102650K .......... .......... .......... .......... .......... 76%  135M 0s
    ## 102700K .......... .......... .......... .......... .......... 76%  125M 0s
    ## 102750K .......... .......... .......... .......... .......... 76%  131M 0s
    ## 102800K .......... .......... .......... .......... .......... 76% 32.9M 0s
    ## 102850K .......... .......... .......... .......... .......... 76% 29.6M 0s
    ## 102900K .......... .......... .......... .......... .......... 76% 91.0M 0s
    ## 102950K .......... .......... .......... .......... .......... 76% 93.5M 0s
    ## 103000K .......... .......... .......... .......... .......... 76%  129M 0s
    ## 103050K .......... .......... .......... .......... .......... 76% 97.2M 0s
    ## 103100K .......... .......... .......... .......... .......... 76% 47.1M 0s
    ## 103150K .......... .......... .......... .......... .......... 76%  115M 0s
    ## 103200K .......... .......... .......... .......... .......... 76%  117M 0s
    ## 103250K .......... .......... .......... .......... .......... 76% 52.9M 0s
    ## 103300K .......... .......... .......... .......... .......... 76% 66.4M 0s
    ## 103350K .......... .......... .......... .......... .......... 76% 49.4M 0s
    ## 103400K .......... .......... .......... .......... .......... 76% 78.2M 0s
    ## 103450K .......... .......... .......... .......... .......... 76%  113M 0s
    ## 103500K .......... .......... .......... .......... .......... 76% 86.4M 0s
    ## 103550K .......... .......... .......... .......... .......... 76% 43.1M 0s
    ## 103600K .......... .......... .......... .......... .......... 76%  105M 0s
    ## 103650K .......... .......... .......... .......... .......... 76%  116M 0s
    ## 103700K .......... .......... .......... .......... .......... 77% 51.6M 0s
    ## 103750K .......... .......... .......... .......... .......... 77% 98.7M 0s
    ## 103800K .......... .......... .......... .......... .......... 77% 93.1M 0s
    ## 103850K .......... .......... .......... .......... .......... 77%  113M 0s
    ## 103900K .......... .......... .......... .......... .......... 77%  123M 0s
    ## 103950K .......... .......... .......... .......... .......... 77%  139M 0s
    ## 104000K .......... .......... .......... .......... .......... 77% 83.7M 0s
    ## 104050K .......... .......... .......... .......... .......... 77% 44.5M 0s
    ## 104100K .......... .......... .......... .......... .......... 77% 72.1M 0s
    ## 104150K .......... .......... .......... .......... .......... 77% 94.2M 0s
    ## 104200K .......... .......... .......... .......... .......... 77% 58.4M 0s
    ## 104250K .......... .......... .......... .......... .......... 77% 52.8M 0s
    ## 104300K .......... .......... .......... .......... .......... 77%  112M 0s
    ## 104350K .......... .......... .......... .......... .......... 77% 33.4M 0s
    ## 104400K .......... .......... .......... .......... .......... 77%  105M 0s
    ## 104450K .......... .......... .......... .......... .......... 77%  126M 0s
    ## 104500K .......... .......... .......... .......... .......... 77%  128M 0s
    ## 104550K .......... .......... .......... .......... .......... 77%  121M 0s
    ## 104600K .......... .......... .......... .......... .......... 77%  123M 0s
    ## 104650K .......... .......... .......... .......... .......... 77% 95.5M 0s
    ## 104700K .......... .......... .......... .......... .......... 77% 91.8M 0s
    ## 104750K .......... .......... .......... .......... .......... 77%  126M 0s
    ## 104800K .......... .......... .......... .......... .......... 77% 41.4M 0s
    ## 104850K .......... .......... .......... .......... .......... 77%  116M 0s
    ## 104900K .......... .......... .......... .......... .......... 77%  116M 0s
    ## 104950K .......... .......... .......... .......... .......... 77% 94.2M 0s
    ## 105000K .......... .......... .......... .......... .......... 77%  103M 0s
    ## 105050K .......... .......... .......... .......... .......... 78% 75.7M 0s
    ## 105100K .......... .......... .......... .......... .......... 78%  127M 0s
    ## 105150K .......... .......... .......... .......... .......... 78% 80.2M 0s
    ## 105200K .......... .......... .......... .......... .......... 78% 39.2M 0s
    ## 105250K .......... .......... .......... .......... .......... 78% 68.3M 0s
    ## 105300K .......... .......... .......... .......... .......... 78%  128M 0s
    ## 105350K .......... .......... .......... .......... .......... 78%  124M 0s
    ## 105400K .......... .......... .......... .......... .......... 78%  125M 0s
    ## 105450K .......... .......... .......... .......... .......... 78%  129M 0s
    ## 105500K .......... .......... .......... .......... .......... 78% 82.7M 0s
    ## 105550K .......... .......... .......... .......... .......... 78% 57.4M 0s
    ## 105600K .......... .......... .......... .......... .......... 78% 67.9M 0s
    ## 105650K .......... .......... .......... .......... .......... 78% 33.5M 0s
    ## 105700K .......... .......... .......... .......... .......... 78% 49.0M 0s
    ## 105750K .......... .......... .......... .......... .......... 78%  130M 0s
    ## 105800K .......... .......... .......... .......... .......... 78%  126M 0s
    ## 105850K .......... .......... .......... .......... .......... 78%  129M 0s
    ## 105900K .......... .......... .......... .......... .......... 78%  132M 0s
    ## 105950K .......... .......... .......... .......... .......... 78% 88.0M 0s
    ## 106000K .......... .......... .......... .......... .......... 78% 67.4M 0s
    ## 106050K .......... .......... .......... .......... .......... 78% 45.9M 0s
    ## 106100K .......... .......... .......... .......... .......... 78% 41.5M 0s
    ## 106150K .......... .......... .......... .......... .......... 78% 47.4M 0s
    ## 106200K .......... .......... .......... .......... .......... 78%  131M 0s
    ## 106250K .......... .......... .......... .......... .......... 78%  129M 0s
    ## 106300K .......... .......... .......... .......... .......... 78%  140M 0s
    ## 106350K .......... .......... .......... .......... .......... 78%  126M 0s
    ## 106400K .......... .......... .......... .......... .......... 79%  110M 0s
    ## 106450K .......... .......... .......... .......... .......... 79% 99.6M 0s
    ## 106500K .......... .......... .......... .......... .......... 79% 82.0M 0s
    ## 106550K .......... .......... .......... .......... .......... 79% 51.4M 0s
    ## 106600K .......... .......... .......... .......... .......... 79% 51.2M 0s
    ## 106650K .......... .......... .......... .......... .......... 79%  113M 0s
    ## 106700K .......... .......... .......... .......... .......... 79%  104M 0s
    ## 106750K .......... .......... .......... .......... .......... 79%  120M 0s
    ## 106800K .......... .......... .......... .......... .......... 79%  114M 0s
    ## 106850K .......... .......... .......... .......... .......... 79%  117M 0s
    ## 106900K .......... .......... .......... .......... .......... 79%  114M 0s
    ## 106950K .......... .......... .......... .......... .......... 79%  121M 0s
    ## 107000K .......... .......... .......... .......... .......... 79%  107M 0s
    ## 107050K .......... .......... .......... .......... .......... 79% 85.3M 0s
    ## 107100K .......... .......... .......... .......... .......... 79%  131M 0s
    ## 107150K .......... .......... .......... .......... .......... 79%  111M 0s
    ## 107200K .......... .......... .......... .......... .......... 79% 70.9M 0s
    ## 107250K .......... .......... .......... .......... .......... 79% 64.3M 0s
    ## 107300K .......... .......... .......... .......... .......... 79% 72.5M 0s
    ## 107350K .......... .......... .......... .......... .......... 79% 85.4M 0s
    ## 107400K .......... .......... .......... .......... .......... 79% 83.8M 0s
    ## 107450K .......... .......... .......... .......... .......... 79% 54.5M 0s
    ## 107500K .......... .......... .......... .......... .......... 79%  113M 0s
    ## 107550K .......... .......... .......... .......... .......... 79% 73.6M 0s
    ## 107600K .......... .......... .......... .......... .......... 79%  119M 0s
    ## 107650K .......... .......... .......... .......... .......... 79% 44.2M 0s
    ## 107700K .......... .......... .......... .......... .......... 79%  104M 0s
    ## 107750K .......... .......... .......... .......... .......... 80%  120M 0s
    ## 107800K .......... .......... .......... .......... .......... 80%  129M 0s
    ## 107850K .......... .......... .......... .......... .......... 80%  127M 0s
    ## 107900K .......... .......... .......... .......... .......... 80%  105M 0s
    ## 107950K .......... .......... .......... .......... .......... 80% 52.9M 0s
    ## 108000K .......... .......... .......... .......... .......... 80% 36.9M 0s
    ## 108050K .......... .......... .......... .......... .......... 80%  126M 0s
    ## 108100K .......... .......... .......... .......... .......... 80% 33.3M 0s
    ## 108150K .......... .......... .......... .......... .......... 80%  106M 0s
    ## 108200K .......... .......... .......... .......... .......... 80%  132M 0s
    ## 108250K .......... .......... .......... .......... .......... 80% 80.3M 0s
    ## 108300K .......... .......... .......... .......... .......... 80% 92.8M 0s
    ## 108350K .......... .......... .......... .......... .......... 80%  131M 0s
    ## 108400K .......... .......... .......... .......... .......... 80%  143M 0s
    ## 108450K .......... .......... .......... .......... .......... 80% 37.8M 0s
    ## 108500K .......... .......... .......... .......... .......... 80%  113M 0s
    ## 108550K .......... .......... .......... .......... .......... 80% 41.9M 0s
    ## 108600K .......... .......... .......... .......... .......... 80% 68.3M 0s
    ## 108650K .......... .......... .......... .......... .......... 80%  128M 0s
    ## 108700K .......... .......... .......... .......... .......... 80% 29.6M 0s
    ## 108750K .......... .......... .......... .......... .......... 80%  149M 0s
    ## 108800K .......... .......... .......... .......... .......... 80%  136M 0s
    ## 108850K .......... .......... .......... .......... .......... 80%  104M 0s
    ## 108900K .......... .......... .......... .......... .......... 80%  127M 0s
    ## 108950K .......... .......... .......... .......... .......... 80%  129M 0s
    ## 109000K .......... .......... .......... .......... .......... 80%  147M 0s
    ## 109050K .......... .......... .......... .......... .......... 80% 46.6M 0s
    ## 109100K .......... .......... .......... .......... .......... 81%  115M 0s
    ## 109150K .......... .......... .......... .......... .......... 81% 61.8M 0s
    ## 109200K .......... .......... .......... .......... .......... 81% 33.4M 0s
    ## 109250K .......... .......... .......... .......... .......... 81% 81.1M 0s
    ## 109300K .......... .......... .......... .......... .......... 81%  127M 0s
    ## 109350K .......... .......... .......... .......... .......... 81%  151M 0s
    ## 109400K .......... .......... .......... .......... .......... 81% 54.1M 0s
    ## 109450K .......... .......... .......... .......... .......... 81%  148M 0s
    ## 109500K .......... .......... .......... .......... .......... 81%  109M 0s
    ## 109550K .......... .......... .......... .......... .......... 81%  103M 0s
    ## 109600K .......... .......... .......... .......... .......... 81%  111M 0s
    ## 109650K .......... .......... .......... .......... .......... 81% 52.9M 0s
    ## 109700K .......... .......... .......... .......... .......... 81% 67.5M 0s
    ## 109750K .......... .......... .......... .......... .......... 81% 52.0M 0s
    ## 109800K .......... .......... .......... .......... .......... 81% 47.4M 0s
    ## 109850K .......... .......... .......... .......... .......... 81%  128M 0s
    ## 109900K .......... .......... .......... .......... .......... 81%  108M 0s
    ## 109950K .......... .......... .......... .......... .......... 81%  118M 0s
    ## 110000K .......... .......... .......... .......... .......... 81%  118M 0s
    ## 110050K .......... .......... .......... .......... .......... 81%  153M 0s
    ## 110100K .......... .......... .......... .......... .......... 81%  130M 0s
    ## 110150K .......... .......... .......... .......... .......... 81%  134M 0s
    ## 110200K .......... .......... .......... .......... .......... 81% 73.6M 0s
    ## 110250K .......... .......... .......... .......... .......... 81%  118M 0s
    ## 110300K .......... .......... .......... .......... .......... 81% 93.9M 0s
    ## 110350K .......... .......... .......... .......... .......... 81%  138M 0s
    ## 110400K .......... .......... .......... .......... .......... 81%  134M 0s
    ## 110450K .......... .......... .......... .......... .......... 82% 79.5M 0s
    ## 110500K .......... .......... .......... .......... .......... 82%  142M 0s
    ## 110550K .......... .......... .......... .......... .......... 82% 75.2M 0s
    ## 110600K .......... .......... .......... .......... .......... 82%  100M 0s
    ## 110650K .......... .......... .......... .......... .......... 82%  126M 0s
    ## 110700K .......... .......... .......... .......... .......... 82% 93.0M 0s
    ## 110750K .......... .......... .......... .......... .......... 82%  166M 0s
    ## 110800K .......... .......... .......... .......... .......... 82% 59.4M 0s
    ## 110850K .......... .......... .......... .......... .......... 82%  111M 0s
    ## 110900K .......... .......... .......... .......... .......... 82% 39.7M 0s
    ## 110950K .......... .......... .......... .......... .......... 82%  126M 0s
    ## 111000K .......... .......... .......... .......... .......... 82%  143M 0s
    ## 111050K .......... .......... .......... .......... .......... 82%  166M 0s
    ## 111100K .......... .......... .......... .......... .......... 82%  107M 0s
    ## 111150K .......... .......... .......... .......... .......... 82% 82.9M 0s
    ## 111200K .......... .......... .......... .......... .......... 82% 36.8M 0s
    ## 111250K .......... .......... .......... .......... .......... 82% 72.9M 0s
    ## 111300K .......... .......... .......... .......... .......... 82% 35.2M 0s
    ## 111350K .......... .......... .......... .......... .......... 82%  128M 0s
    ## 111400K .......... .......... .......... .......... .......... 82%  109M 0s
    ## 111450K .......... .......... .......... .......... .......... 82%  119M 0s
    ## 111500K .......... .......... .......... .......... .......... 82% 40.1M 0s
    ## 111550K .......... .......... .......... .......... .......... 82% 81.5M 0s
    ## 111600K .......... .......... .......... .......... .......... 82%  123M 0s
    ## 111650K .......... .......... .......... .......... .......... 82%  115M 0s
    ## 111700K .......... .......... .......... .......... .......... 82%  140M 0s
    ## 111750K .......... .......... .......... .......... .......... 82%  133M 0s
    ## 111800K .......... .......... .......... .......... .......... 83% 71.3M 0s
    ## 111850K .......... .......... .......... .......... .......... 83% 54.7M 0s
    ## 111900K .......... .......... .......... .......... .......... 83% 51.3M 0s
    ## 111950K .......... .......... .......... .......... .......... 83%  130M 0s
    ## 112000K .......... .......... .......... .......... .......... 83% 34.6M 0s
    ## 112050K .......... .......... .......... .......... .......... 83%  155M 0s
    ## 112100K .......... .......... .......... .......... .......... 83%  101M 0s
    ## 112150K .......... .......... .......... .......... .......... 83%  147M 0s
    ## 112200K .......... .......... .......... .......... .......... 83%  142M 0s
    ## 112250K .......... .......... .......... .......... .......... 83% 22.5M 0s
    ## 112300K .......... .......... .......... .......... .......... 83%  142M 0s
    ## 112350K .......... .......... .......... .......... .......... 83% 52.7M 0s
    ## 112400K .......... .......... .......... .......... .......... 83%  138M 0s
    ## 112450K .......... .......... .......... .......... .......... 83%  157M 0s
    ## 112500K .......... .......... .......... .......... .......... 83% 68.4M 0s
    ## 112550K .......... .......... .......... .......... .......... 83%  137M 0s
    ## 112600K .......... .......... .......... .......... .......... 83% 58.7M 0s
    ## 112650K .......... .......... .......... .......... .......... 83% 86.6M 0s
    ## 112700K .......... .......... .......... .......... .......... 83% 49.2M 0s
    ## 112750K .......... .......... .......... .......... .......... 83%  115M 0s
    ## 112800K .......... .......... .......... .......... .......... 83% 95.0M 0s
    ## 112850K .......... .......... .......... .......... .......... 83%  180M 0s
    ## 112900K .......... .......... .......... .......... .......... 83% 95.7M 0s
    ## 112950K .......... .......... .......... .......... .......... 83% 85.5M 0s
    ## 113000K .......... .......... .......... .......... .......... 83% 75.1M 0s
    ## 113050K .......... .......... .......... .......... .......... 83% 70.5M 0s
    ## 113100K .......... .......... .......... .......... .......... 83%  137M 0s
    ## 113150K .......... .......... .......... .......... .......... 84% 67.0M 0s
    ## 113200K .......... .......... .......... .......... .......... 84% 83.3M 0s
    ## 113250K .......... .......... .......... .......... .......... 84% 79.4M 0s
    ## 113300K .......... .......... .......... .......... .......... 84% 77.4M 0s
    ## 113350K .......... .......... .......... .......... .......... 84% 81.1M 0s
    ## 113400K .......... .......... .......... .......... .......... 84%  122M 0s
    ## 113450K .......... .......... .......... .......... .......... 84%  103M 0s
    ## 113500K .......... .......... .......... .......... .......... 84% 99.0M 0s
    ## 113550K .......... .......... .......... .......... .......... 84% 23.3M 0s
    ## 113600K .......... .......... .......... .......... .......... 84% 76.2M 0s
    ## 113650K .......... .......... .......... .......... .......... 84% 26.5M 0s
    ## 113700K .......... .......... .......... .......... .......... 84% 70.8M 0s
    ## 113750K .......... .......... .......... .......... .......... 84% 60.3M 0s
    ## 113800K .......... .......... .......... .......... .......... 84% 67.1M 0s
    ## 113850K .......... .......... .......... .......... .......... 84% 87.2M 0s
    ## 113900K .......... .......... .......... .......... .......... 84% 77.1M 0s
    ## 113950K .......... .......... .......... .......... .......... 84% 68.7M 0s
    ## 114000K .......... .......... .......... .......... .......... 84%  117M 0s
    ## 114050K .......... .......... .......... .......... .......... 84%  109M 0s
    ## 114100K .......... .......... .......... .......... .......... 84%  112M 0s
    ## 114150K .......... .......... .......... .......... .......... 84%  108M 0s
    ## 114200K .......... .......... .......... .......... .......... 84% 88.2M 0s
    ## 114250K .......... .......... .......... .......... .......... 84% 84.1M 0s
    ## 114300K .......... .......... .......... .......... .......... 84% 95.4M 0s
    ## 114350K .......... .......... .......... .......... .......... 84% 55.2M 0s
    ## 114400K .......... .......... .......... .......... .......... 84% 76.4M 0s
    ## 114450K .......... .......... .......... .......... .......... 84% 77.7M 0s
    ## 114500K .......... .......... .......... .......... .......... 85% 71.8M 0s
    ## 114550K .......... .......... .......... .......... .......... 85%  131M 0s
    ## 114600K .......... .......... .......... .......... .......... 85% 74.3M 0s
    ## 114650K .......... .......... .......... .......... .......... 85% 79.5M 0s
    ## 114700K .......... .......... .......... .......... .......... 85% 89.5M 0s
    ## 114750K .......... .......... .......... .......... .......... 85% 68.0M 0s
    ## 114800K .......... .......... .......... .......... .......... 85% 88.0M 0s
    ## 114850K .......... .......... .......... .......... .......... 85% 78.3M 0s
    ## 114900K .......... .......... .......... .......... .......... 85% 66.0M 0s
    ## 114950K .......... .......... .......... .......... .......... 85% 80.7M 0s
    ## 115000K .......... .......... .......... .......... .......... 85% 93.1M 0s
    ## 115050K .......... .......... .......... .......... .......... 85% 60.4M 0s
    ## 115100K .......... .......... .......... .......... .......... 85% 82.0M 0s
    ## 115150K .......... .......... .......... .......... .......... 85% 91.2M 0s
    ## 115200K .......... .......... .......... .......... .......... 85% 95.1M 0s
    ## 115250K .......... .......... .......... .......... .......... 85% 71.6M 0s
    ## 115300K .......... .......... .......... .......... .......... 85% 93.8M 0s
    ## 115350K .......... .......... .......... .......... .......... 85% 81.6M 0s
    ## 115400K .......... .......... .......... .......... .......... 85% 86.2M 0s
    ## 115450K .......... .......... .......... .......... .......... 85% 61.3M 0s
    ## 115500K .......... .......... .......... .......... .......... 85%  102M 0s
    ## 115550K .......... .......... .......... .......... .......... 85%  105M 0s
    ## 115600K .......... .......... .......... .......... .......... 85% 83.8M 0s
    ## 115650K .......... .......... .......... .......... .......... 85% 97.1M 0s
    ## 115700K .......... .......... .......... .......... .......... 85% 97.6M 0s
    ## 115750K .......... .......... .......... .......... .......... 85% 65.5M 0s
    ## 115800K .......... .......... .......... .......... .......... 85% 85.9M 0s
    ## 115850K .......... .......... .......... .......... .......... 86% 85.6M 0s
    ## 115900K .......... .......... .......... .......... .......... 86% 90.0M 0s
    ## 115950K .......... .......... .......... .......... .......... 86% 84.4M 0s
    ## 116000K .......... .......... .......... .......... .......... 86%  100M 0s
    ## 116050K .......... .......... .......... .......... .......... 86% 91.9M 0s
    ## 116100K .......... .......... .......... .......... .......... 86% 85.0M 0s
    ## 116150K .......... .......... .......... .......... .......... 86%  115M 0s
    ## 116200K .......... .......... .......... .......... .......... 86%  102M 0s
    ## 116250K .......... .......... .......... .......... .......... 86% 75.3M 0s
    ## 116300K .......... .......... .......... .......... .......... 86% 68.7M 0s
    ## 116350K .......... .......... .......... .......... .......... 86%  117M 0s
    ## 116400K .......... .......... .......... .......... .......... 86% 95.5M 0s
    ## 116450K .......... .......... .......... .......... .......... 86%  106M 0s
    ## 116500K .......... .......... .......... .......... .......... 86% 10.3M 0s
    ## 116550K .......... .......... .......... .......... .......... 86%  126M 0s
    ## 116600K .......... .......... .......... .......... .......... 86%  142M 0s
    ## 116650K .......... .......... .......... .......... .......... 86%  144M 0s
    ## 116700K .......... .......... .......... .......... .......... 86%  131M 0s
    ## 116750K .......... .......... .......... .......... .......... 86%  155M 0s
    ## 116800K .......... .......... .......... .......... .......... 86%  133M 0s
    ## 116850K .......... .......... .......... .......... .......... 86% 60.1M 0s
    ## 116900K .......... .......... .......... .......... .......... 86% 74.2M 0s
    ## 116950K .......... .......... .......... .......... .......... 86% 36.4M 0s
    ## 117000K .......... .......... .......... .......... .......... 86%  111M 0s
    ## 117050K .......... .......... .......... .......... .......... 86%  131M 0s
    ## 117100K .......... .......... .......... .......... .......... 86% 68.0M 0s
    ## 117150K .......... .......... .......... .......... .......... 86%  139M 0s
    ## 117200K .......... .......... .......... .......... .......... 87%  129M 0s
    ## 117250K .......... .......... .......... .......... .......... 87%  147M 0s
    ## 117300K .......... .......... .......... .......... .......... 87% 48.9M 0s
    ## 117350K .......... .......... .......... .......... .......... 87% 77.5M 0s
    ## 117400K .......... .......... .......... .......... .......... 87% 97.6M 0s
    ## 117450K .......... .......... .......... .......... .......... 87% 26.6M 0s
    ## 117500K .......... .......... .......... .......... .......... 87%  130M 0s
    ## 117550K .......... .......... .......... .......... .......... 87% 46.2M 0s
    ## 117600K .......... .......... .......... .......... .......... 87%  127M 0s
    ## 117650K .......... .......... .......... .......... .......... 87%  147M 0s
    ## 117700K .......... .......... .......... .......... .......... 87%  134M 0s
    ## 117750K .......... .......... .......... .......... .......... 87% 18.2M 0s
    ## 117800K .......... .......... .......... .......... .......... 87% 96.7M 0s
    ## 117850K .......... .......... .......... .......... .......... 87%  165M 0s
    ## 117900K .......... .......... .......... .......... .......... 87% 92.2M 0s
    ## 117950K .......... .......... .......... .......... .......... 87%  128M 0s
    ## 118000K .......... .......... .......... .......... .......... 87%  132M 0s
    ## 118050K .......... .......... .......... .......... .......... 87%  150M 0s
    ## 118100K .......... .......... .......... .......... .......... 87%  116M 0s
    ## 118150K .......... .......... .......... .......... .......... 87% 24.2M 0s
    ## 118200K .......... .......... .......... .......... .......... 87% 42.7M 0s
    ## 118250K .......... .......... .......... .......... .......... 87% 77.5M 0s
    ## 118300K .......... .......... .......... .......... .......... 87% 98.5M 0s
    ## 118350K .......... .......... .......... .......... .......... 87% 95.2M 0s
    ## 118400K .......... .......... .......... .......... .......... 87%  142M 0s
    ## 118450K .......... .......... .......... .......... .......... 87% 87.5M 0s
    ## 118500K .......... .......... .......... .......... .......... 87% 72.9M 0s
    ## 118550K .......... .......... .......... .......... .......... 88%  137M 0s
    ## 118600K .......... .......... .......... .......... .......... 88%  101M 0s
    ## 118650K .......... .......... .......... .......... .......... 88%  121M 0s
    ## 118700K .......... .......... .......... .......... .......... 88% 41.9M 0s
    ## 118750K .......... .......... .......... .......... .......... 88%  111M 0s
    ## 118800K .......... .......... .......... .......... .......... 88%  103M 0s
    ## 118850K .......... .......... .......... .......... .......... 88%  118M 0s
    ## 118900K .......... .......... .......... .......... .......... 88%  140M 0s
    ## 118950K .......... .......... .......... .......... .......... 88%  138M 0s
    ## 119000K .......... .......... .......... .......... .......... 88%  108M 0s
    ## 119050K .......... .......... .......... .......... .......... 88%  139M 0s
    ## 119100K .......... .......... .......... .......... .......... 88%  111M 0s
    ## 119150K .......... .......... .......... .......... .......... 88% 64.9M 0s
    ## 119200K .......... .......... .......... .......... .......... 88% 41.8M 0s
    ## 119250K .......... .......... .......... .......... .......... 88%  129M 0s
    ## 119300K .......... .......... .......... .......... .......... 88%  124M 0s
    ## 119350K .......... .......... .......... .......... .......... 88%  119M 0s
    ## 119400K .......... .......... .......... .......... .......... 88%  135M 0s
    ## 119450K .......... .......... .......... .......... .......... 88% 83.4M 0s
    ## 119500K .......... .......... .......... .......... .......... 88%  152M 0s
    ## 119550K .......... .......... .......... .......... .......... 88% 98.1M 0s
    ## 119600K .......... .......... .......... .......... .......... 88% 43.0M 0s
    ## 119650K .......... .......... .......... .......... .......... 88%  139M 0s
    ## 119700K .......... .......... .......... .......... .......... 88% 57.5M 0s
    ## 119750K .......... .......... .......... .......... .......... 88% 89.3M 0s
    ## 119800K .......... .......... .......... .......... .......... 88% 99.9M 0s
    ## 119850K .......... .......... .......... .......... .......... 88% 43.5M 0s
    ## 119900K .......... .......... .......... .......... .......... 89% 65.7M 0s
    ## 119950K .......... .......... .......... .......... .......... 89%  164M 0s
    ## 120000K .......... .......... .......... .......... .......... 89%  132M 0s
    ## 120050K .......... .......... .......... .......... .......... 89%  153M 0s
    ## 120100K .......... .......... .......... .......... .......... 89%  101M 0s
    ## 120150K .......... .......... .......... .......... .......... 89%  107M 0s
    ## 120200K .......... .......... .......... .......... .......... 89%  120M 0s
    ## 120250K .......... .......... .......... .......... .......... 89%  150M 0s
    ## 120300K .......... .......... .......... .......... .......... 89% 60.7M 0s
    ## 120350K .......... .......... .......... .......... .......... 89%  148M 0s
    ## 120400K .......... .......... .......... .......... .......... 89%  134M 0s
    ## 120450K .......... .......... .......... .......... .......... 89%  111M 0s
    ## 120500K .......... .......... .......... .......... .......... 89%  121M 0s
    ## 120550K .......... .......... .......... .......... .......... 89% 45.4M 0s
    ## 120600K .......... .......... .......... .......... .......... 89%  114M 0s
    ## 120650K .......... .......... .......... .......... .......... 89%  159M 0s
    ## 120700K .......... .......... .......... .......... .......... 89% 67.5M 0s
    ## 120750K .......... .......... .......... .......... .......... 89% 51.1M 0s
    ## 120800K .......... .......... .......... .......... .......... 89% 75.4M 0s
    ## 120850K .......... .......... .......... .......... .......... 89%  124M 0s
    ## 120900K .......... .......... .......... .......... .......... 89% 24.6M 0s
    ## 120950K .......... .......... .......... .......... .......... 89%  153M 0s
    ## 121000K .......... .......... .......... .......... .......... 89% 54.3M 0s
    ## 121050K .......... .......... .......... .......... .......... 89%  139M 0s
    ## 121100K .......... .......... .......... .......... .......... 89% 99.0M 0s
    ## 121150K .......... .......... .......... .......... .......... 89%  135M 0s
    ## 121200K .......... .......... .......... .......... .......... 89%  131M 0s
    ## 121250K .......... .......... .......... .......... .......... 90%  155M 0s
    ## 121300K .......... .......... .......... .......... .......... 90% 85.2M 0s
    ## 121350K .......... .......... .......... .......... .......... 90%  115M 0s
    ## 121400K .......... .......... .......... .......... .......... 90%  115M 0s
    ## 121450K .......... .......... .......... .......... .......... 90%  165M 0s
    ## 121500K .......... .......... .......... .......... .......... 90%  128M 0s
    ## 121550K .......... .......... .......... .......... .......... 90%  142M 0s
    ## 121600K .......... .......... .......... .......... .......... 90% 58.0M 0s
    ## 121650K .......... .......... .......... .......... .......... 90% 89.6M 0s
    ## 121700K .......... .......... .......... .......... .......... 90% 63.9M 0s
    ## 121750K .......... .......... .......... .......... .......... 90% 54.4M 0s
    ## 121800K .......... .......... .......... .......... .......... 90% 39.8M 0s
    ## 121850K .......... .......... .......... .......... .......... 90% 63.2M 0s
    ## 121900K .......... .......... .......... .......... .......... 90%  105M 0s
    ## 121950K .......... .......... .......... .......... .......... 90% 65.9M 0s
    ## 122000K .......... .......... .......... .......... .......... 90% 25.1M 0s
    ## 122050K .......... .......... .......... .......... .......... 90%  145M 0s
    ## 122100K .......... .......... .......... .......... .......... 90%  121M 0s
    ## 122150K .......... .......... .......... .......... .......... 90%  130M 0s
    ## 122200K .......... .......... .......... .......... .......... 90%  133M 0s
    ## 122250K .......... .......... .......... .......... .......... 90%  156M 0s
    ## 122300K .......... .......... .......... .......... .......... 90%  104M 0s
    ## 122350K .......... .......... .......... .......... .......... 90%  155M 0s
    ## 122400K .......... .......... .......... .......... .......... 90% 40.5M 0s
    ## 122450K .......... .......... .......... .......... .......... 90% 70.3M 0s
    ## 122500K .......... .......... .......... .......... .......... 90% 27.2M 0s
    ## 122550K .......... .......... .......... .......... .......... 90%  108M 0s
    ## 122600K .......... .......... .......... .......... .......... 91%  124M 0s
    ## 122650K .......... .......... .......... .......... .......... 91%  112M 0s
    ## 122700K .......... .......... .......... .......... .......... 91% 86.5M 0s
    ## 122750K .......... .......... .......... .......... .......... 91%  113M 0s
    ## 122800K .......... .......... .......... .......... .......... 91% 99.8M 0s
    ## 122850K .......... .......... .......... .......... .......... 91% 56.5M 0s
    ## 122900K .......... .......... .......... .......... .......... 91% 99.2M 0s
    ## 122950K .......... .......... .......... .......... .......... 91% 95.7M 0s
    ## 123000K .......... .......... .......... .......... .......... 91% 49.3M 0s
    ## 123050K .......... .......... .......... .......... .......... 91% 62.9M 0s
    ## 123100K .......... .......... .......... .......... .......... 91%  113M 0s
    ## 123150K .......... .......... .......... .......... .......... 91%  122M 0s
    ## 123200K .......... .......... .......... .......... .......... 91% 82.7M 0s
    ## 123250K .......... .......... .......... .......... .......... 91%  125M 0s
    ## 123300K .......... .......... .......... .......... .......... 91% 79.5M 0s
    ## 123350K .......... .......... .......... .......... .......... 91% 80.9M 0s
    ## 123400K .......... .......... .......... .......... .......... 91% 48.2M 0s
    ## 123450K .......... .......... .......... .......... .......... 91%  108M 0s
    ## 123500K .......... .......... .......... .......... .......... 91%  114M 0s
    ## 123550K .......... .......... .......... .......... .......... 91%  102M 0s
    ## 123600K .......... .......... .......... .......... .......... 91% 51.1M 0s
    ## 123650K .......... .......... .......... .......... .......... 91% 87.5M 0s
    ## 123700K .......... .......... .......... .......... .......... 91%  117M 0s
    ## 123750K .......... .......... .......... .......... .......... 91%  119M 0s
    ## 123800K .......... .......... .......... .......... .......... 91% 60.5M 0s
    ## 123850K .......... .......... .......... .......... .......... 91% 86.2M 0s
    ## 123900K .......... .......... .......... .......... .......... 91% 99.4M 0s
    ## 123950K .......... .......... .......... .......... .......... 92%  116M 0s
    ## 124000K .......... .......... .......... .......... .......... 92% 50.1M 0s
    ## 124050K .......... .......... .......... .......... .......... 92%  111M 0s
    ## 124100K .......... .......... .......... .......... .......... 92% 88.0M 0s
    ## 124150K .......... .......... .......... .......... .......... 92%  136M 0s
    ## 124200K .......... .......... .......... .......... .......... 92% 94.4M 0s
    ## 124250K .......... .......... .......... .......... .......... 92% 79.4M 0s
    ## 124300K .......... .......... .......... .......... .......... 92%  111M 0s
    ## 124350K .......... .......... .......... .......... .......... 92% 73.0M 0s
    ## 124400K .......... .......... .......... .......... .......... 92%  116M 0s
    ## 124450K .......... .......... .......... .......... .......... 92% 70.2M 0s
    ## 124500K .......... .......... .......... .......... .......... 92% 74.2M 0s
    ## 124550K .......... .......... .......... .......... .......... 92%  109M 0s
    ## 124600K .......... .......... .......... .......... .......... 92%  109M 0s
    ## 124650K .......... .......... .......... .......... .......... 92% 57.5M 0s
    ## 124700K .......... .......... .......... .......... .......... 92%  118M 0s
    ## 124750K .......... .......... .......... .......... .......... 92%  102M 0s
    ## 124800K .......... .......... .......... .......... .......... 92% 98.0M 0s
    ## 124850K .......... .......... .......... .......... .......... 92% 67.4M 0s
    ## 124900K .......... .......... .......... .......... .......... 92% 99.7M 0s
    ## 124950K .......... .......... .......... .......... .......... 92%  119M 0s
    ## 125000K .......... .......... .......... .......... .......... 92% 90.2M 0s
    ## 125050K .......... .......... .......... .......... .......... 92% 83.5M 0s
    ## 125100K .......... .......... .......... .......... .......... 92%  115M 0s
    ## 125150K .......... .......... .......... .......... .......... 92% 76.5M 0s
    ## 125200K .......... .......... .......... .......... .......... 92% 87.7M 0s
    ## 125250K .......... .......... .......... .......... .......... 92% 84.3M 0s
    ## 125300K .......... .......... .......... .......... .......... 93%  105M 0s
    ## 125350K .......... .......... .......... .......... .......... 93%  112M 0s
    ## 125400K .......... .......... .......... .......... .......... 93%  109M 0s
    ## 125450K .......... .......... .......... .......... .......... 93% 87.1M 0s
    ## 125500K .......... .......... .......... .......... .......... 93% 97.6M 0s
    ## 125550K .......... .......... .......... .......... .......... 93% 94.2M 0s
    ## 125600K .......... .......... .......... .......... .......... 93%  122M 0s
    ## 125650K .......... .......... .......... .......... .......... 93% 66.3M 0s
    ## 125700K .......... .......... .......... .......... .......... 93%  100M 0s
    ## 125750K .......... .......... .......... .......... .......... 93%  131M 0s
    ## 125800K .......... .......... .......... .......... .......... 93% 64.8M 0s
    ## 125850K .......... .......... .......... .......... .......... 93% 98.3M 0s
    ## 125900K .......... .......... .......... .......... .......... 93%  102M 0s
    ## 125950K .......... .......... .......... .......... .......... 93% 92.1M 0s
    ## 126000K .......... .......... .......... .......... .......... 93%  105M 0s
    ## 126050K .......... .......... .......... .......... .......... 93% 99.4M 0s
    ## 126100K .......... .......... .......... .......... .......... 93%  116M 0s
    ## 126150K .......... .......... .......... .......... .......... 93% 70.4M 0s
    ## 126200K .......... .......... .......... .......... .......... 93% 95.8M 0s
    ## 126250K .......... .......... .......... .......... .......... 93% 70.1M 0s
    ## 126300K .......... .......... .......... .......... .......... 93% 80.3M 0s
    ## 126350K .......... .......... .......... .......... .......... 93%  132M 0s
    ## 126400K .......... .......... .......... .......... .......... 93% 56.6M 0s
    ## 126450K .......... .......... .......... .......... .......... 93%  124M 0s
    ## 126500K .......... .......... .......... .......... .......... 93%  119M 0s
    ## 126550K .......... .......... .......... .......... .......... 93%  162M 0s
    ## 126600K .......... .......... .......... .......... .......... 93%  121M 0s
    ## 126650K .......... .......... .......... .......... .......... 94% 26.4M 0s
    ## 126700K .......... .......... .......... .......... .......... 94%  126M 0s
    ## 126750K .......... .......... .......... .......... .......... 94% 50.9M 0s
    ## 126800K .......... .......... .......... .......... .......... 94% 44.6M 0s
    ## 126850K .......... .......... .......... .......... .......... 94% 74.0M 0s
    ## 126900K .......... .......... .......... .......... .......... 94%  104M 0s
    ## 126950K .......... .......... .......... .......... .......... 94%  123M 0s
    ## 127000K .......... .......... .......... .......... .......... 94%  142M 0s
    ## 127050K .......... .......... .......... .......... .......... 94% 25.4M 0s
    ## 127100K .......... .......... .......... .......... .......... 94%  157M 0s
    ## 127150K .......... .......... .......... .......... .......... 94%  170M 0s
    ## 127200K .......... .......... .......... .......... .......... 94%  148M 0s
    ## 127250K .......... .......... .......... .......... .......... 94%  167M 0s
    ## 127300K .......... .......... .......... .......... .......... 94%  149M 0s
    ## 127350K .......... .......... .......... .......... .......... 94%  124M 0s
    ## 127400K .......... .......... .......... .......... .......... 94% 66.7M 0s
    ## 127450K .......... .......... .......... .......... .......... 94% 60.1M 0s
    ## 127500K .......... .......... .......... .......... .......... 94% 50.8M 0s
    ## 127550K .......... .......... .......... .......... .......... 94%  151M 0s
    ## 127600K .......... .......... .......... .......... .......... 94%  155M 0s
    ## 127650K .......... .......... .......... .......... .......... 94%  107M 0s
    ## 127700K .......... .......... .......... .......... .......... 94%  108M 0s
    ## 127750K .......... .......... .......... .......... .......... 94%  137M 0s
    ## 127800K .......... .......... .......... .......... .......... 94% 31.1M 0s
    ## 127850K .......... .......... .......... .......... .......... 94%  113M 0s
    ## 127900K .......... .......... .......... .......... .......... 94% 39.3M 0s
    ## 127950K .......... .......... .......... .......... .......... 94%  137M 0s
    ## 128000K .......... .......... .......... .......... .......... 95% 96.4M 0s
    ## 128050K .......... .......... .......... .......... .......... 95% 39.7M 0s
    ## 128100K .......... .......... .......... .......... .......... 95% 95.7M 0s
    ## 128150K .......... .......... .......... .......... .......... 95% 84.7M 0s
    ## 128200K .......... .......... .......... .......... .......... 95% 67.5M 0s
    ## 128250K .......... .......... .......... .......... .......... 95%  109M 0s
    ## 128300K .......... .......... .......... .......... .......... 95%  113M 0s
    ## 128350K .......... .......... .......... .......... .......... 95%  169M 0s
    ## 128400K .......... .......... .......... .......... .......... 95%  142M 0s
    ## 128450K .......... .......... .......... .......... .......... 95%  130M 0s
    ## 128500K .......... .......... .......... .......... .......... 95%  126M 0s
    ## 128550K .......... .......... .......... .......... .......... 95% 94.7M 0s
    ## 128600K .......... .......... .......... .......... .......... 95% 54.9M 0s
    ## 128650K .......... .......... .......... .......... .......... 95% 66.2M 0s
    ## 128700K .......... .......... .......... .......... .......... 95% 41.9M 0s
    ## 128750K .......... .......... .......... .......... .......... 95% 61.0M 0s
    ## 128800K .......... .......... .......... .......... .......... 95% 75.1M 0s
    ## 128850K .......... .......... .......... .......... .......... 95%  161M 0s
    ## 128900K .......... .......... .......... .......... .......... 95% 65.2M 0s
    ## 128950K .......... .......... .......... .......... .......... 95%  153M 0s
    ## 129000K .......... .......... .......... .......... .......... 95%  132M 0s
    ## 129050K .......... .......... .......... .......... .......... 95%  155M 0s
    ## 129100K .......... .......... .......... .......... .......... 95% 60.3M 0s
    ## 129150K .......... .......... .......... .......... .......... 95% 69.9M 0s
    ## 129200K .......... .......... .......... .......... .......... 95% 49.0M 0s
    ## 129250K .......... .......... .......... .......... .......... 95%  129M 0s
    ## 129300K .......... .......... .......... .......... .......... 95% 38.7M 0s
    ## 129350K .......... .......... .......... .......... .......... 96%  140M 0s
    ## 129400K .......... .......... .......... .......... .......... 96% 51.4M 0s
    ## 129450K .......... .......... .......... .......... .......... 96%  164M 0s
    ## 129500K .......... .......... .......... .......... .......... 96%  128M 0s
    ## 129550K .......... .......... .......... .......... .......... 96%  181M 0s
    ## 129600K .......... .......... .......... .......... .......... 96% 59.2M 0s
    ## 129650K .......... .......... .......... .......... .......... 96%  121M 0s
    ## 129700K .......... .......... .......... .......... .......... 96% 74.0M 0s
    ## 129750K .......... .......... .......... .......... .......... 96% 86.4M 0s
    ## 129800K .......... .......... .......... .......... .......... 96% 25.5M 0s
    ## 129850K .......... .......... .......... .......... .......... 96% 48.9M 0s
    ## 129900K .......... .......... .......... .......... .......... 96% 80.5M 0s
    ## 129950K .......... .......... .......... .......... .......... 96%  127M 0s
    ## 130000K .......... .......... .......... .......... .......... 96%  129M 0s
    ## 130050K .......... .......... .......... .......... .......... 96%  168M 0s
    ## 130100K .......... .......... .......... .......... .......... 96%  105M 0s
    ## 130150K .......... .......... .......... .......... .......... 96%  162M 0s
    ## 130200K .......... .......... .......... .......... .......... 96%  128M 0s
    ## 130250K .......... .......... .......... .......... .......... 96% 57.5M 0s
    ## 130300K .......... .......... .......... .......... .......... 96% 33.6M 0s
    ## 130350K .......... .......... .......... .......... .......... 96% 68.2M 0s
    ## 130400K .......... .......... .......... .......... .......... 96% 61.1M 0s
    ## 130450K .......... .......... .......... .......... .......... 96%  132M 0s
    ## 130500K .......... .......... .......... .......... .......... 96% 42.3M 0s
    ## 130550K .......... .......... .......... .......... .......... 96% 74.6M 0s
    ## 130600K .......... .......... .......... .......... .......... 96%  129M 0s
    ## 130650K .......... .......... .......... .......... .......... 97%  128M 0s
    ## 130700K .......... .......... .......... .......... .......... 97% 57.5M 0s
    ## 130750K .......... .......... .......... .......... .......... 97% 66.5M 0s
    ## 130800K .......... .......... .......... .......... .......... 97%  118M 0s
    ## 130850K .......... .......... .......... .......... .......... 97% 69.1M 0s
    ## 130900K .......... .......... .......... .......... .......... 97% 38.6M 0s
    ## 130950K .......... .......... .......... .......... .......... 97%  123M 0s
    ## 131000K .......... .......... .......... .......... .......... 97% 97.4M 0s
    ## 131050K .......... .......... .......... .......... .......... 97%  111M 0s
    ## 131100K .......... .......... .......... .......... .......... 97% 41.0M 0s
    ## 131150K .......... .......... .......... .......... .......... 97%  107M 0s
    ## 131200K .......... .......... .......... .......... .......... 97% 78.4M 0s
    ## 131250K .......... .......... .......... .......... .......... 97% 49.3M 0s
    ## 131300K .......... .......... .......... .......... .......... 97% 68.2M 0s
    ## 131350K .......... .......... .......... .......... .......... 97%  153M 0s
    ## 131400K .......... .......... .......... .......... .......... 97% 85.6M 0s
    ## 131450K .......... .......... .......... .......... .......... 97% 80.6M 0s
    ## 131500K .......... .......... .......... .......... .......... 97% 64.6M 0s
    ## 131550K .......... .......... .......... .......... .......... 97%  127M 0s
    ## 131600K .......... .......... .......... .......... .......... 97% 70.8M 0s
    ## 131650K .......... .......... .......... .......... .......... 97% 52.1M 0s
    ## 131700K .......... .......... .......... .......... .......... 97%  101M 0s
    ## 131750K .......... .......... .......... .......... .......... 97%  113M 0s
    ## 131800K .......... .......... .......... .......... .......... 97% 74.3M 0s
    ## 131850K .......... .......... .......... .......... .......... 97% 51.5M 0s
    ## 131900K .......... .......... .......... .......... .......... 97% 94.4M 0s
    ## 131950K .......... .......... .......... .......... .......... 97%  157M 0s
    ## 132000K .......... .......... .......... .......... .......... 98% 64.4M 0s
    ## 132050K .......... .......... .......... .......... .......... 98% 56.0M 0s
    ## 132100K .......... .......... .......... .......... .......... 98%  115M 0s
    ## 132150K .......... .......... .......... .......... .......... 98% 72.2M 0s
    ## 132200K .......... .......... .......... .......... .......... 98%  116M 0s
    ## 132250K .......... .......... .......... .......... .......... 98% 56.7M 0s
    ## 132300K .......... .......... .......... .......... .......... 98% 95.2M 0s
    ## 132350K .......... .......... .......... .......... .......... 98%  128M 0s
    ## 132400K .......... .......... .......... .......... .......... 98% 67.8M 0s
    ## 132450K .......... .......... .......... .......... .......... 98% 67.7M 0s
    ## 132500K .......... .......... .......... .......... .......... 98%  125M 0s
    ## 132550K .......... .......... .......... .......... .......... 98% 99.1M 0s
    ## 132600K .......... .......... .......... .......... .......... 98% 49.6M 0s
    ## 132650K .......... .......... .......... .......... .......... 98% 90.0M 0s
    ## 132700K .......... .......... .......... .......... .......... 98% 89.3M 0s
    ## 132750K .......... .......... .......... .......... .......... 98% 97.5M 0s
    ## 132800K .......... .......... .......... .......... .......... 98% 46.3M 0s
    ## 132850K .......... .......... .......... .......... .......... 98%  101M 0s
    ## 132900K .......... .......... .......... .......... .......... 98%  111M 0s
    ## 132950K .......... .......... .......... .......... .......... 98%  126M 0s
    ## 133000K .......... .......... .......... .......... .......... 98% 78.3M 0s
    ## 133050K .......... .......... .......... .......... .......... 98% 81.5M 0s
    ## 133100K .......... .......... .......... .......... .......... 98%  103M 0s
    ## 133150K .......... .......... .......... .......... .......... 98%  121M 0s
    ## 133200K .......... .......... .......... .......... .......... 98% 46.9M 0s
    ## 133250K .......... .......... .......... .......... .......... 98% 89.4M 0s
    ## 133300K .......... .......... .......... .......... .......... 98% 93.5M 0s
    ## 133350K .......... .......... .......... .......... .......... 99%  150M 0s
    ## 133400K .......... .......... .......... .......... .......... 99%  117M 0s
    ## 133450K .......... .......... .......... .......... .......... 99% 71.7M 0s
    ## 133500K .......... .......... .......... .......... .......... 99% 89.5M 0s
    ## 133550K .......... .......... .......... .......... .......... 99% 92.8M 0s
    ## 133600K .......... .......... .......... .......... .......... 99% 66.3M 0s
    ## 133650K .......... .......... .......... .......... .......... 99%  146M 0s
    ## 133700K .......... .......... .......... .......... .......... 99% 54.3M 0s
    ## 133750K .......... .......... .......... .......... .......... 99%  128M 0s
    ## 133800K .......... .......... .......... .......... .......... 99% 70.4M 0s
    ## 133850K .......... .......... .......... .......... .......... 99% 73.4M 0s
    ## 133900K .......... .......... .......... .......... .......... 99%  121M 0s
    ## 133950K .......... .......... .......... .......... .......... 99%  125M 0s
    ## 134000K .......... .......... .......... .......... .......... 99% 72.0M 0s
    ## 134050K .......... .......... .......... .......... .......... 99%  143M 0s
    ## 134100K .......... .......... .......... .......... .......... 99% 65.9M 0s
    ## 134150K .......... .......... .......... .......... .......... 99%  107M 0s
    ## 134200K .......... .......... .......... .......... .......... 99% 92.7M 0s
    ## 134250K .......... .......... .......... .......... .......... 99% 83.6M 0s
    ## 134300K .......... .......... .......... .......... .......... 99% 86.1M 0s
    ## 134350K .......... .......... .......... .......... .......... 99%  131M 0s
    ## 134400K .......... .......... .......... .......... .......... 99% 88.0M 0s
    ## 134450K .......... .......... .......... .......... .......... 99% 80.4M 0s
    ## 134500K .......... .......... .......... .......... .......... 99% 89.0M 0s
    ## 134550K .......... .......... .......... .......... .......... 99%  146M 0s
    ## 134600K .......... .......... .......... .......... .......... 99% 70.6M 0s
    ## 134650K .......... .......... .......... .......... .......... 99%  111M 0s
    ## 134700K .......... .......... .......... ..........           100%  128M=1.9s
    ## 
    ## 2020-12-19 14:08:44 (68.2 MB/s) - ‘silva_nr99_v138_train_set.fa.gz.3’ saved [137973851/137973851]

``` r
taxa <- assignTaxonomy(seqtab.nochim, "~/CC2_Dada2AndPhyloseq/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
```

``` r
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
```

    ##      Kingdom    Phylum             Class                 Order            
    ## [1,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [2,] "Bacteria" "Cyanobacteria"    "Cyanobacteriia"      "Synechococcales"
    ## [3,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [4,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [5,] "Bacteria" "Proteobacteria"   "Alphaproteobacteria" "SAR11 clade"    
    ## [6,] "Bacteria" "Actinobacteriota" "Acidimicrobiia"      "Actinomarinales"
    ##      Family             Genus                    
    ## [1,] "Clade I"          "Clade Ia"               
    ## [2,] "Cyanobiaceae"     "Synechococcus CC9902"   
    ## [3,] "Clade I"          "Clade Ia"               
    ## [4,] "Clade I"          "Clade Ia"               
    ## [5,] "Clade II"         NA                       
    ## [6,] "Actinomarinaceae" "Candidatus Actinomarina"

# Phyloseq

## Importation des librairies et préparation des données

``` r
library(phyloseq)
library(Biostrings)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:phyloseq':
    ## 
    ##     distance

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(ggplot2)
```

``` r
samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "Y"), `[`, 1)
echantillon <- substr(subject,10,999)
subject <- substr(subject,1,999)
date <- sapply(strsplit(samples.out, "_"), `[`, 3)
profondeur <- sapply(strsplit(samples.out, "_"), `[`, 2)
samdf <- data.frame(Subject=subject, Echantillon=echantillon, Date=date, Profondeur=profondeur)
samdf$date <- "Sep14" [samdf$Date!="Sep14"] 
samdf$date <- "Mars15" [samdf$Date!="Mars15"] 
rownames(samdf) <- samples.out
```

``` r
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
```

Construction de l’objet PhyloSeq : on associe la table des OTU, les
données des échantillons et de la table taxonomique à ps.

``` r
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 1557 taxa and 11 samples ]
    ## sample_data() Sample Data:       [ 11 samples by 5 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 1557 taxa by 6 taxonomic ranks ]
    ## refseq()      DNAStringSet:      [ 1557 reference sequences ]

## Visualisation de l’alpha diversité

``` r
plot_richness(ps, x="Profondeur", measures=c("Shannon", "Simpson"), color="Date")
```

![](02_dada2_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

La date du 10 Septembre sera considérée comme celle de l’été tandis que
le 11 Mars 2015, comme celle de l’hiver. A partir de ces graphes, on
peut voir qu’il n’y a pas de grande différence entre les communautés
microbiennes des fonds marins entre les deux dates. Pour les
échantillons Median, il n’y a pas d’éléments à la date du 11 Mars 2015
mais on peut voir que les deux échatillons sont proches entre eux.
Enfin, les échantillons de Surface sont nettement différents entre les
deux périodes. Cela peut s’expliquer pour l’intensité lumineuse, qui est
plus forte et plus longue en été qu’en hiver et que certaines bactéries
(comme les cyanoblaceae que nous verrons plus tard) sont bien plus
présentes. De même pour la température, corrélée à l’intensité
lumineuse, qui est plus forte en été et dont seules les bactéries
mésophiles pourront croître. Concernant les fonds marins, on n’observe
pas de grandes différences car les bactéries présentes sont très peu
affectées par les changements de saisons. On retrouvera les mêmes
bactéries présentes continuellement au cours de l’année, celles ne
nécessitant pas d’énergie lumineuse, d’oxygène ou pouvant croître à des
températures assez faibles.

## Ordination

``` r
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray")
```

``` r
plot_ordination(ps.prop, ord.pcoa.bray, shape = "Date", color="Profondeur", title="PCoA Bray")
```

![](02_dada2_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

Question 1 : On retrouve les communautés de Surface aux niveaux des
extrêmités gauches et droites. On retrouve les communautés Median
isolées à gauche, au niveau de -0.2. Enfin, pour les communautés Fond,
on les retrouves en haut du graphe et vers 0.3 de l’axe 1. D’après les
résultats du graphique, on observe une corrélation concernant certains
échantillons du 10 Septembre 2014 qui forment trois groupes bien
distinsts, tandis que les échantillons du 11 Mars 2015 sont isolés entre
eux. D’après ces premiers résultats, on peut confirmer que la saison a
un impact important sur la composition de la communauté microbienne.
Concernant sur la profondeur, on observe que les échantillons d’hiver
(11 Mars 2015) ne sont pas très différents entre eux, ils sont mêmes
similaires, à quelques différences près. Par contre, pour les
échantillons d’été (10 Septembre 2014), les différentes sont bien plus
nettes. Le groupe isolé en haut, représentant les fonds marins, est bien
différent des échantillons de surface et median. Ainsi, on peut
confirmer les hypothèses proposées dans la partie sur les indices de
similarité sont vraies.

## Représention des échantillons en “histogramme” (bar plot)

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Profondeur", fill="Family") + facet_wrap(~Date, scales="free_x")
```

![](02_dada2_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

Question 2 : On peut voir, qu’entre ces deux graphes, la famille des
Cyanoblaceae est présente en 2014 mais pas, ou très peu, en 2015. De
même pour la famille des Rhodobacteraceaea. Ces deux familles peuvent
constituer des biomarqueurs pour savoir comment la communauté
microbienne évolue et déduire quelles bactéries seront plus présentes
que d’autres en fonction des saisons.

``` r
top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Profondeur", fill="Genus") + facet_wrap(~Date, scales="free_x")
```

![](02_dada2_files/figure-gfm/unnamed-chunk-33-1.png)<!-- -->

En regardant les genres présents dans les familles, on peut voir que ce
sont les Synechococcus CC9902, de la famille des Cyanoblaceae, qui
prédominent en été mais qui sont quasi absent en hiver.
