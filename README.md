**Ioan** is a reverse annotator for resolving ambiguous mutation annotations.

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

#### Infer codon identity
Given two amino acid positions and infer potential identity due to different usage of transcripts.

```
#!bash
ioan codonanno -a hg19.map -c 
```

 + input: 1) codon position 1; 2) codon position 2;



## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.