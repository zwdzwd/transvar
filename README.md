**Ioan** is a reverse annotator for resolving ambiguous mutation annotations. Ioan can read files in their .gz format.

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

#### Find nucleotide position(s) given amino acid positions

```
#!bash
ioan codonanno -a hg19.map -c PIK3CA:E545K
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
$ ioan codoneq -c MET.p1010 MET.p992 --ensembl ~/reference/ensembl/Homo_sapiens.GRCh37.75.gtf.gz --ref ~/reference/hs37d5.fa

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