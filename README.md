**Ioan** is a reverse annotator for inferring genomic characterization(s) of mutational events (e.g., ```chr3:178936091 G=>A```) from transcript-dependent mutation annotation(s) (e.g., ```PIK3CA p.E545K```). It is designed for resolving ambiguous mutation annotations arising from differential transcript usage. Ioan supports transcript annotation from commonly-used databases including Ensembl, NCBI RefSeq, GENCODE, CCDS, UCSC etc. Ioan can read files in their .gz format.

--------

[TOC]

--------

### Download and Install

#### program
```
#!bash

 $ wget https://bitbucket.org/wanding/ioan/get/v1.0.zip
 $ unzip [downloaded zip]
 $ make
```

#### transcript annotations

##### Ensembl

in GTF format, can be downloaded from 

http://http://www.ensembl.org/info/data/ftp/index.html

can be used in ioan via ```--ensembl Homo_sapiens.GRCh37.75.gtf.gz```

##### NCBI RefSeq

in GFF3 format, can be downloaded from

ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz

can be used in ioan via option ```--refseq ref_GRCh37.p13_top_level.gff3.gz```

##### CCDS



### Usage


#### reverse annotate amino acid mutation

```
#!bash
ioan codonanno --ref ~/reference/hs37d5.fa
    --ccds ~/reference/CCDS/Hs37.3/CCDS.current.txt
    -c PIK3CA:E545K
```
outputs
```
PIK3CA:E545K    CDDS    CCDS43171.1     PIK3CA  545 \
 3   178936091   178936092   178936093   GAG  +   E=>K \
178936091  G  A       AAG,AAA
```
 + input: 1) transcript annotation file; 2) codon position; 3) (optional) mutation information;
 + output: 1) annotation;


#### reverse annotate nucleotide mutation

```
#!bash
ioan nucanno -c BRAF.c1781A>G
```


#### infer potential codon identity
Given two amino acid positions and infer potential identity due to different usage of transcripts.

```
#!bash
$ ioan codoneq -c MET.p1010 MET.p992
     --ensembl ~/reference/ensembl/Homo_sapiens.GRCh37.75.gtf.gz
     --ref ~/reference/hs37d5.fa
```
gives
```
MET 1010
transcript [ENST00000397752] 1  codon: 116412043-116414935-116414936
transcript [ENST00000318493] 2  codon: 116411989-116411990-116411991
MET 992
transcript [ENST00000397752] 1  codon: 116411989,116411990,116411991
transcript [ENST00000318493] 2  codon: 116411935,116411936,116411937
Genomic location might be the same.

```

 + input: 1) gene-codon position 1; 2) gene-codon position 2; 3) annotation database

#### search alternative codon identifiers
Given a codon identifier, search the transcript annotations for alternative (codon) identifiers
```
ioan codonsearch 
```

#### annotate genomic mutations

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.