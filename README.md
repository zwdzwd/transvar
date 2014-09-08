**Ioan** is a reverse annotator for inferring genomic characterization(s) of mutations (e.g., ```chr3:178936091 G=>A```) from transcript-dependent annotation(s) (e.g., ```PIK3CA p.E545K``` or ```PIK3CA c.1633G>A```, which are extensively used in clinical settings). It is designed for resolving ambiguous mutation annotations arising from differential transcript usage. Ioan supports transcript annotation from commonly-used databases including Ensembl, NCBI RefSeq, GENCODE, CCDS, UCSC etc. Ioan can read files in their .gz format. Ioan provides functionality of forward annotation as well.

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

The following list the transcript annotations supported by ioan. Ioan can take any one or any combination(s) of these transcript annotations as long as these annotations are based on the same version of reference assembly.

##### Ensembl (GTF)

 + available [here](http://http://www.ensembl.org/info/data/ftp/index.html)

 + used in ioan via option ```--ensembl Homo_sapiens.GRCh37.75.gtf.gz```

##### NCBI RefSeq (GFF3)

 + available [here](ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz)

 + used in ioan via option ```--refseq ref_GRCh37.p13_top_level.gff3.gz```

##### CCDS (table)

 + available [here](http://www.ncbi.nlm.nih.gov/CCDS/CcdsBrowse.cgi)

 + used in ioan via option ```--ccds CCDS.current.txt```

##### GENCODE (GTF)

 + available [here](http://www.gencodegenes.org/releases/19.html)

 + used in ioan via option ```--gencode gencode.v19.annotation.gtf.gz```

##### UCSC knownGene (table)

 + available [here](https://genome.ucsc.edu/cgi-bin/hgTables?command=start)

 + used in ioan via option ```--kg UCSC_knownGene_hg19.gz --alias UCSC_kgAlias.gz```

##### RefGene via UCSC (table)

 + available [here](https://genome.ucsc.edu/cgi-bin/hgTables?command=start)

 + used in ioan via option ```--ucsc2 hg19.map```

### Usage


#### reverse annotation of amino acid mutations
Ioan automatically recognizes the amino acid mutations. Acceptable mutation formats are ```PIK3CA:E545K``` or ```PIK3CA:p.E545K```, or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.

```
#!bash
$ ioan revanno --ref hs37d5.fa --ccds CCDS.current.txt -i PIK3CA:E545K
```
outputs
```
#!text
PIK3CA:E545K       CDDS       CCDS43171.1       PIK3CA 
    545    3     178936091     178936092     178936093
    GAG    +     E=>K    178936091    G    A       AAG,AAA
```
 + input: 1) transcript annotation file; 2) codon position; 3) (optional) mutation information;
 + output: 1) annotation;


#### reverse annotation of nucleotide mutations
Ioan infers nucleotide mutation through ```PIK3CA:1633G>A``` or ```PIK3CA:c.1633G>A```. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.
```
#!bash
$ ioan revanno --ref hs37d5.fa --ccds CCDS.current.txt -i 'PIK3CA:c.1633G>A'
```
outputs

```
#!text
PIK3CA:c.1633G>A        CDDS    CCDS43171.1     PIK3CA
    545     3       178936091       178936092       178936093
    GAG     +       E=>K    178936091       G       A
```

#### infer potential codon identity
Given two amino acid positions and infer potential identity due to different usage of transcripts.

```
#!bash
$ ioan codoneq -c MET.p1010 MET.p992 --ensembl Homo_sapiens.GRCh37.75.gtf.gz --ref hs37d5.fa
```
gives
```
#!text
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
#!bash
$ ioan codonsearch -l input.table --ccds CCDS.current.txt --ref hs37d5.fa -m 1 -o 1
```
Example input.table
```
#!text
CDKN2A.p58
CDKN2A.p61
CDKN2A.p69
CDKN2A.p69
ERBB2.p755
ERBB2.p755
```
will output
```
#!text
CDKN2A.p58      73      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p61      76      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      84      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      55      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
ERBB2.p755      785     CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
ERBB2.p755      725     CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
```

#### annotate mutations with genomic locations (forward annotation)

```
#!bash
ioan anno --ccds CCDS.current.txt --ref hs37d5.fa -i 'chr3:178936091.G>A'
```
outputs
```
#!text
chr3:178936091.G>A      CDDS    CCDS43171.1     PIK3CA  545     3       178936091       178936092       178936093       GAG     +       E=>K    178936091       G     A
```

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.