**Ioan** is a reverse annotator for resolving ambiguous mutation annotations. It is for discovering discrepancy in annotation due to difference in transcript usage. Ioan can read files in their .gz format.

--------

[TOC]

--------

### Download and Install

```
#!bash

 $ wget https://bitbucket.org/wanding/ioan/get/v1.0.zip
 $ unzip [downloaded zip]
 $ make
```

### Usage

#### search alternative identifiers


#### search a list of mutations accounting of potential difference in transcript usage


#### Find nucleotide position(s) given amino acid positions

```
#!bash
ioan codonanno --ref ~/reference/hs37d5.fa
    --ccds ~/reference/CCDS/Hs37.3/CCDS.current.txt
    -c PIK3CA:E545K
```
outputs
```
PIK3CA:E545K    CDDS    CCDS43171.1     PIK3CA  545 \
 3   178936091   178936092   178936093   GAG  +   E=>K 178936091  G  A       AAG,AAA
```
 + input: 1) transcript annotation file; 2) codon position; 3) (optional) mutation information;
 + output: 1) annotation;

#### Find amino acid position(s) given amino acid positions

```
#!bash
ioan codonsearch -a hg19.map -c PIK3CA:E545K
```

#### Infer potential codon identity
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

#### Search list of mutations (codon level) for matching target mutation (codon level)
This corresponds to, for example, when one needs to search a database of recurrent mutations or driver/hotspot mutations for a target mutation.


```
#!bash
ioan codonmatch -a hg19.map -c
```


## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.