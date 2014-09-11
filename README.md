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
Ioan automatically recognizes the amino acid mutations. Acceptable mutation formats are ```PIK3CA:E545K``` or ```PIK3CA:p.E545K```, or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. Ioan takes native HGVS format inputs and outputs. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.

```
#!bash
$ ioan revanno -i PIK3CA:E545K --ref hs37d5.fa --ensembl Homo_sapiens.GRCh37.75.gtf.gz
```
outputs
```
#!text
PIK3CA:E545K    3       178936091-178936092-178936093   ENST00000263967 
PIK3CA (+, coding)      3:G178936091A/c.1633G>A/p.E545K 
CddMuts=3:G178936091A;NCodonSeq=GAG;NCddSeqs=AAG,AAA
```
 + input: 1) transcript annotation file; 2) codon position; 3) (optional) mutation information;
 + output: 1) annotation;

RevAn may encounter cases where the ambiguity cannot be fully resolved. For example,
```
#!text
ACSL4   23:108926078    108926078       c.399C>T        p.R133R Missense
      X       108926078-108926079-108926080   CCDS14548.1      ACSL4 (-, coding)       X:G108926078A/c.399C>T/p.R133R
      CddMuts=X:G108926078T,X:G108926078C,X:G108926078A;NCodonSeq=CGC;NCddSeqs=AGG,AGA,CGA,CGC,CGG,CGT
```
In those cases, RevAn prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

#### reverse annotation of nucleotide mutations
Ioan infers nucleotide mutation through ```PIK3CA:1633G>A``` or ```PIK3CA:c.1633G>A```. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.
```
#!bash
$ ioan revanno --ref hs37d5.fa --ccds CCDS.current.txt -i 'PIK3CA:c.1633G>A'
```
outputs

```
#!text
PIK3CA:E545K    3       178936091-178936092-178936093   CCDS43171.1
     PIK3CA (+, coding)      3:G178936091A/c.1633G>A/p.E545K
 CddMuts=3:G178936091A;NCodonSeq=GAG;NCddSeqs=AAG,AAA
```
or one can batch process a list of mutation identifiers with optional transcript id to constraint the search
```
#!bash
$ ioan revanno -l input_table -g 1 -m 4 -t 2 --ensembl Homo_sapiens.GRCh37.75.gtf.gz --ref hs37d5.fa -o 2,3,5 | les
```
As suggested by the command, RevAn takes as input the 1st column as gene and 4th column as identifier. The 2nd column will be used as the transcript id from Ensembl to constrain the alternative identifier search. The 2nd, 3rd and 5th columns are chosen to be output as a validation of RevAn's performance.
```
#!text
ADAMTSL3        ENST00000286744 15:84442328     c.243G>A        p.W81*  Nonsense
ADAMTSL3        ENST00000286744 15:84442326     c.241T>C        p.W81R  Missense
ADAMTSL4        ENST00000369038 1:150530513     c.2270G>A       p.G757D Missense
ADCY2   ENST00000338316 5:7802364       c.2662G>A       p.V888I Missense
ADCY2   ENST00000338316 5:7802365       c.2663T>C       p.V888A Missense
```
and output
```
#!text
ENST00000286744 15:84442328     p.W81*  15      84442326-84442327-84442328      ENST00000286744 ADAMTSL3 (+ coding)     15:G84442328A/c.243G>A/p.W81*   CddMuts=15:G84442328A;NCodonSeq=TGG;NCddSeqs=TGA
ENST00000286744 15:84442326     p.W81R  15      84442326-84442327-84442328      ENST00000286744 ADAMTSL3 (+ coding)     15:T84442326C/c.241T>C/p.W81R   CddMuts=15:T84442326C;NCodonSeq=TGG;NCddSeqs=CGG
ENST00000369038 1:150530513     p.G757D 1       150530512-150530513-150530514   ENST00000369038 ADAMTSL4 (+ coding)     1:G150530513A/c.2270G>A/p.G757D CddMuts=1:G150530513A;NCodonSeq=GGT;NCddSeqs=GAT
ENST00000338316 5:7802364       p.V888I 5       7802364-7802365-7802366 ENST00000338316 ADCY2 (+ coding)        5:G7802364A/c.2662G>A/p.V888I   CddMuts=5:G7802364A;NCodonSeq=GTC;NCddSeqs=ATC
ENST00000338316 5:7802365       p.V888A 5       7802364-7802365-7802366 ENST00000338316 ADCY2 (+ coding)        5:T7802365C/c.2663T>C/p.V888A   CddMuts=5:T7802365C;NCodonSeq=GTC;NCddSeqs=GCC
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
$ ioan codonsearch -i DHODH:G152R --ref ~/reference/hs37d5.fa --refseq ~/reference/refseq/ref_GRCh37.p13_top_level.gff3.gz
```
outputs
```
#!text
DHODH:G152R     p.G9    16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq]
DHODH:G152R     p.G124  16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq]
DHODH:G152R     p.G16   16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq]
```
RevAn outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

One can also run `ioan codonsearch` to batch process a list of mutation identifiers.
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
outputs
```
#!text
CDKN2A.p58      p.72    9       21971184-21971185-21971186      21971185-21971186-21971187      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p58      p.73    9       21971184-21971185-21971186      21971182-21971183-21971184      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p61      p.75    9       21971175-21971176-21971177      21971176-21971177-21971178      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p61      p.76    9       21971175-21971176-21971177      21971173-21971174-21971175      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      p.55    9       21971194-21971195-21971196      21971193-21971194-21971195      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS],CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      p.84    9       21971151-21971152-21971153      21971149-21971150-21971151      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      p.83    9       21971151-21971152-21971153      21971152-21971153-21971154      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      p.54    9       21971194-21971195-21971196      21971196-21971197-21971198      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      p.55    9       21971194-21971195-21971196      21971193-21971194-21971195      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS],CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      p.84    9       21971151-21971152-21971153      21971149-21971150-21971151      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      p.83    9       21971151-21971152-21971153      21971152-21971153-21971154      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      p.54    9       21971194-21971195-21971196      21971196-21971197-21971198      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
ERBB2.p755      p.725   17      37880219-37880220-37880221      37880219-37880220-37880221      CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
ERBB2.p755      p.785   17      37881024-37881025-37881026      37881024-37881025-37881026      CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
ERBB2.p755      p.725   17      37880219-37880220-37880221      37880219-37880220-37880221      CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
ERBB2.p755      p.785   17      37881024-37881025-37881026      37881024-37881025-37881026      CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
```
The third column indicates the potential transcript usage for the alternative identifier. Each transcript usage is denoted by <listing transcript>/<actual transcript>. Different potential choices are separated by ','.

#### annotate mutations from genomic locations

This is the forward annotation

```
#!bash
ioan anno --ccds CCDS.current.txt --ref hs37d5.fa -i 'chr3:178936091.G>A'
```
outputs
```
#!text
chr3:178936091.G>A      3       178936091-178936092-178936093   CCDS43171.1
     PIK3CA (+, coding)      3:G178936091A/c.1633>/p.E545K   .
```

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.