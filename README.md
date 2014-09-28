**RevAn** is a reverse annotator for inferring genomic characterization(s) of mutations (e.g., ```chr3:178936091 G=>A```) from transcript-dependent annotation(s) (e.g., ```PIK3CA p.E545K``` or ```PIK3CA c.1633G>A```). It is designed for resolving ambiguous mutation annotations arising from differential transcript usage. RevAn has the following features:

 + supports uncertainty in HGVS nomenclature
 + supports single nucleotide variation (SNV), insertions and deletions (indels) and block substitutions
 + supports mutations at both coding region and intronic/UTR regions
 + supports transcript annotation from commonly-used databases such as Ensembl, NCBI RefSeq and GENCODE etc
 + functionality of forward annotation.

--------

[TOC]

--------

### Quick start

```
#!bash
 $ wget https://bitbucket.org/wanding/revan/get/v1.01.zip
 $ unzip v1.01.zip
 $ cd [unzipped dir]
 $ ./revan config --download_hg19_anno
 $ ./revan revanno --ucsc -i 'PIK3CA.p.E545K'
```

### Download and Install

#### program
```
#!bash

 $ wget https://bitbucket.org/wanding/revan/get/v1.0.zip
 $ unzip [downloaded zip]
 $ make
```

#### reference genome assembly
For most annotation database (the only exception may be the UCSC table which is used as example in the Quick Start), RevAn requires a samtools indexed reference genome in fasta format, which is available at, e.g., [UCSC ftp](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/).
Once downloaded and indexed, one could use the reference in RevAn through the `--ref` option followed by the fasta filename. One also has the option of specifying the default location to `revan.cfg` by:
```
#!bash
revan config -k reference -v hg19.fa
```
so that there is no need to specify the location of reference on subsequent usages.

#### transcript annotations

The following list the transcript annotations supported by revan. RevAn can take any one or any combination(s) of these transcript annotations as long as these annotations are based on the same version of reference assembly.

```
#!bash
revan config --download_hg19_anno
```
will automatically download annotation from
Ensembl, RefSeq, UCSC RefGene, GENCODE, AceView and UCSC knownGene to `[install dir]/download` directory.

### Usage

#### specify transcript annotation

The following table summarize the option(s) to use each database in the annotation

 | Database | Format | Default  | Non-default     |
 |:-------|:--------|:---------|:----|
 | CCDS | CCDS flat text | `--ccds` | `--ccds CCDS.current.txt` |
 | UCSC RefGene | UCSC flat text | `--ucsc` | `--ucsc2 hg19.ucsc.txt` |
 | Ensembl | Ensembl GTF | `--ensembl`  | `--ensembl GRCh37.75.gtf.gz`  |
 | RefSeq | RefSeq GFF3 | `--refseq`  | `--refseq GRCh37.p13.gff3.gz`   |
 | AceView | AceView GFF | `--aceview` | `--aceview AceView.ncbi37.gff.gz`  |
 | GENCODE | GENCODE GTF | `--gencode` | `--gencode gencode.v19.gtf.gz`  |
 | knownGene | knownGene table | `-kg` | `--kg kg.gz --alias kgAlias.gz` |

If one download transcripts through `revan config`, RevAn would use the downloaded definition automatically (by setting the default configuration file). For example, `--ccds` would look for the downloaded CCDS definition. One can specify non-default annotation by `--ccds [annotation file]`. Or one can set the default annotation by 
```revan config -k ccds -v [annotation file]```. 
The configuration file is located either at the `[install dir]/revan.cfg` or `~/.revan.cfg` if the installation directory is inaccessible.

---

#### batch processing

For all mutation types, one can batch process a list of mutation identifiers with optional transcript id to constraint the search. Take SNV for example,
```
#!bash
$ revan revanno -l input_table -g 1 -m 4 -t 2 --ensembl -o 2,3,5
```
As suggested by the command, RevAn takes as input the 1st column as gene and 4th column as identifier. The 2nd column will be used as the transcript id from Ensembl to constrain the alternative identifier search. The 2nd, 3rd and 5th columns are chosen to be output as a validation of RevAn's performance.

Input:
```
#!text
ADAMTSL3        ENST00000286744 15:84442328     c.243G>A        p.W81*  Nonsense
ADAMTSL3        ENST00000286744 15:84442326     c.241T>C        p.W81R  Missense
ADAMTSL4        ENST00000369038 1:150530513     c.2270G>A       p.G757D Missense
ADCY2   ENST00000338316 5:7802364       c.2662G>A       p.V888I Missense
ADCY2   ENST00000338316 5:7802365       c.2663T>C       p.V888A Missense
```
Output:
```
#!text
ENST00000286744 15:84442328     p.W81*  15      84442326-84442327-84442328      ENST00000286744 ADAMTSL3 (+ coding)     15:G84442328A/c.243G>A/p.W81*   CddMuts=15:G84442328A;NCodonSeq=TGG;NCddSeqs=TGA
ENST00000286744 15:84442326     p.W81R  15      84442326-84442327-84442328      ENST00000286744 ADAMTSL3 (+ coding)     15:T84442326C/c.241T>C/p.W81R   CddMuts=15:T84442326C;NCodonSeq=TGG;NCddSeqs=CGG
ENST00000369038 1:150530513     p.G757D 1       150530512-150530513-150530514   ENST00000369038 ADAMTSL4 (+ coding)     1:G150530513A/c.2270G>A/p.G757D CddMuts=1:G150530513A;NCodonSeq=GGT;NCddSeqs=GAT
ENST00000338316 5:7802364       p.V888I 5       7802364-7802365-7802366 ENST00000338316 ADCY2 (+ coding)        5:G7802364A/c.2662G>A/p.V888I   CddMuts=5:G7802364A;NCodonSeq=GTC;NCddSeqs=ATC
ENST00000338316 5:7802365       p.V888A 5       7802364-7802365-7802366 ENST00000338316 ADCY2 (+ coding)        5:T7802365C/c.2663T>C/p.V888A   CddMuts=5:T7802365C;NCodonSeq=GTC;NCddSeqs=GCC
```

---

#### reverse annotation of single amino acid substitution
Mutation formats acceptable in RevAn are ```PIK3CA:E545K``` or ```PIK3CA:p.E545K```, or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. RevAn takes native HGVS format inputs and outputs. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.

```
#!bash
$ revan revanno -i PIK3CA:E545K --ensembl
```
outputs
```
#!text
PIK3CA:E545K    3       178936091-178936092-178936093   ENST00000263967
    PIK3CA (+, coding)      3:G178936091A/c.1633G>A/p.E545K 
    CddMuts=3:G178936091A;NCodonSeq=GAG;NCddSeqs=AAG,AAA
```

One may encounter **ambiguous cases** where the multiple substitutions exist in explaining the amino acid change. For example,
```
#!text
ACSL4   23:108926078    108926078       c.399C>T        p.R133R Missense
    X       108926078-108926079-108926080   CCDS14548.1
    ACSL4 (-, coding)       X:G108926078A/c.399C>T/p.R133R
    CddMuts=X:G108926078T,X:G108926078C,X:G108926078A;NCodonSeq=CGC;NCddSeqs=AGG,AGA,CGA,CGC,CGG,CGT
```
In those cases, RevAn prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

---
#### reverse annotation of single nucleotide variation (SNV)

RevAn infers nucleotide mutation through ```PIK3CA:1633G>A``` or ```PIK3CA:c.1633G>A```. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.
```
#!bash
$ revan revanno --ccds -i 'PIK3CA:c.1633G>A'
```
outputs
```
#!text
PIK3CA:E545K    3    178936091-178936092-178936093   CCDS43171.1
    PIK3CA (+, coding)      3:G178936091A/c.1633G>A/p.E545K
    CddMuts=3:G178936091A;NCodonSeq=GAG;NCddSeqs=AAG,AAA
```
---

#### reverse annotation of nucleotide insertion
An insertion may result in: 1) a pure insertion of amino acids; 2) a block substitution of amino acids, when insertion occur after 1st or 2nd base in a codon; or 3) a frame-shift. Following HGVS nomenclature, RevAn labels the first different amino acid and the length of the peptide util stop codon, assuming no change in the splicing.

Example: to annotate a **in-frame, in-phase insertion**,
```
#!bash
$ revan revanno --ccds -i 'ACIN1:c.1932_1933insATTCAC'
```
```
#!text
ACIN1:c.1932_1933insATTCAC    14     14:23548785-(ins)-23548786
    CCDS55905.1    ACIN1 (-, coding)       
    14:23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH
    NatInsSeq=ATTCAC;RefInsSeq=GTGAAT;Phase=0
```
`Phase = 0,1,2` indicates whether the insertion happen after the 3rd, 1st or 2nd base of a codon, respectively. An insertion *in phase* refers to one with `Phase=0`.

Example: to annotate an **out-of-phase, in-frame insertion**,
```
#!bash
$ revan revanno --ccds -i 'ACIN1:c.1930_1931insATTCAC'
```
```
#!text
ACIN1:c.1930_1931insATTCAC   14     14:23548787-(ins)-23548788      CCDS9587.1
    ACIN1 (-, coding)       14:23548787_23548788insGTGAAT/c.1930_1931insATTCAC/p.S643_R644insHS
    NatInsSeq=C(ATTCAC)GT;RefInsSeq=GTGAAT;Phase=1
```

Example: to annotate a **frame-shift insertion**,
```
#!bash
$ revan revanno --ccds -i 'AAAS:c.1225_1226insG'
```
```
#!text
AAAS:c.1225_1226insG    12   12:53702089-53702090    CCDS8856.1    AAAS (-, coding)
    12:53702089_53702090insC/c.1225_1226insG/p.E409Gfs*17   NatInsSeq=G;RefInsSeq=C
```

Example: to annotate an **intronic insertion**,
```
#!bash
$ revan revanno --ccds -i 'ADAM33:c.991-3_991-2insC'
```
```
#!text
ADAM33:c.991-3_991-2insC   20   20:3654141-3654142-3654143-(3654145)-(ins)-(3654146)
    CCDS13058.1     ADAM33 (- intronic)     20:3654145_3654146insG/c.991-3_991-2insC/.
    RefInsSeq=G;NatInsSeq=C
```
In the case of intronic insertions, amino acid identifier is not applicable, represented in a `.`.

---

#### reverse annotation of nucleotide deletion
Similar to insertions, deletion can be in-frame or frame-shift. The consequence of deletion to amino acid sequence may appear a simple deletion or a block substitution (in the case where in-frame deletion is out of phase, i.e., partially delete codons).

Example: to annotate an **in-frame deletion**,
```
#!bash
$ revan revanno --ccds -i 'A4GNT:c.694_696delTTG'
```
```
#!text
A4GNT:c.694_696delTTG   3   137843433-137843435     CCDS3097.1    A4GNT (- coding)
    3:137843433_137843435del/c.694_696del/p.L232del RefDelSeq=CAA;NatDelSeq=TTG
```

Example: to annotate a **in-frame, out-of-phase deletion**,
```
#!bash
$ revan revanno --ccds -i 'ABHD15.c.431_433delGTG'
```
```
#!text
ABHD15.c.431_433delGTG  17   27893552-27893554    CCDS32602.1     ABHD15 (- coding)
    17:27893552_27893554del/c.431_433del/p.C144_V145delinsF RefDelSeq=CAC;NatDelSeq=GTG
```

Example: to annotate a **frame-shift deletion**,
```
#!bash
$ revan revanno --ccds -i 'AADACL3.c.374delG'
```
```
#!text
AADACL3.c.374delG    1    12785494-12785494    CCDS41252.1   AADACL3 (+ coding)
    1:12785494_12785494del/c.374del/p.C125Ffs*17    RefDelSeq=G;NatDelSeq=G
```

Example: to annotate a **deletion that span from intronic to coding region**,
```
#!bash
$ revan revanno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'
```
```
#!text
ABCB11:c.1198-8_1199delcactccagAA       2       2:169833196-169833205
   CCDS46444.1     ABCB11 (- coding & intronic)
    2:169833196_169833205del/c.1198-8_1199del/p.K400Tfs*4
   RefDelSeq=TTCTGGAGTG;NatDelSeq=CACTCCAGAA
```

---

#### reverse annotation of nucleotide block substitution

Example: to annotate a block substitution in **coding region**,
```
#!bash
$ revan revanno --ccds -i 'A1CF:c.508_509CC>TT'
```
```
#!text
A1CF:c.508_509CC>TT     10      52595929-52595930       CCDS7241.1
      A1CF (-, coding)        10:52595929_52595930GG>AA/c.508_509CC>TT/p.P170L        .
A1CF:c.508_509CC>TT     10      52595929-52595930       CCDS7242.1
      A1CF (-, coding)        10:52595929_52595930GG>AA/c.508_509CC>TT/p.P170L        .
```

Block substitution does not necessarily results in block substitution in amino acid. For example, the following substitution results in a deletion.
```
#!bash
$ revan revanno --ccds -i 'CSRNP1.c.1212_1224>GGAGGAGGAA'
```
```
#!text
CSRNP1.c.1212_1224>GGAGGAGGAA   3    39185092-39185104   CCDS2682.1   CSRNP1 (-, coding)
    3:39185092_39185104TTCCTCCTCCTCC>TTCCTCCTCC/c.1212_1224GGAGGAGGAGGAA>GGAGGAGGAA/p.E408del .
```

Example: to annotate a block substitution in **intronic region**,
```
#!bash
$ revan revanno --ccds -i 'A1CF:c.1460+2_1460+3TG>CC'
```
```
#!text
A1CF:c.1460+2_1460+3TG>CC    10    52570797-52570798   CCDS7241.1    A1CF (-, intronic)
    10:52570797_52570798CA>GG/c.1460+2_1460+3TG>CC/.        .
```

When block substitution occurs **across splice site**, RevAn put a tag in the info fields and does not predict amino acid change.
```
#!bash
$ revan revanno --ccds -i 'A1CF:c.1459_1460+3ATGTG>CC'
```
```
#!text
A1CF:c.1459_1460+3ATGTG>CC    10   52570797-52570801   CCDS7241.1   A1CF (-, coding;intronic)
    10:52570797_52570801CACAT>GG/c.1459_1460+3ATGTG>CC/.    CrossSplitSite
```

---
#### infer potential codon identity
Given two amino acid positions and infer potential identity due to different usage of transcripts.

```
#!bash
$ revan codoneq -c MET.p1010 MET.p992 --ensembl Homo_sapiens.GRCh37.75.gtf.gz --ref hs37d5.fa
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

Example: to search alternative identifiers
```
#!bash
$ revan codonsearch --ccds -i CDKN2A.p.58
```
```
#!text
CDKN2A.p.58    CDKN2A.p.72   9    21971184-21971185-21971186   21971185-21971186-21971187
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p.58    CDKN2A.p.73   9    21971184-21971185-21971186   21971182-21971183-21971184
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
```

Example2: to search alternative identifiers of DHODH:G152R
```
#!bash
$ revan codonsearch -i DHODH:G152R --refseq
```
outputs
```
#!text
DHODH:G152R     DHODH.p.G9    16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq]
DHODH:G152R     DHODH.p.G124  16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq]
DHODH:G152R     DHODH.p.G16   16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq]
```
RevAn outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

Example: to run `revan codonsearch` to batch process a list of mutation identifiers.
```
#!bash
$ revan codonsearch -l input.table --ccds -m 1 -o 1
```
Example input.table
```
#!text
CDKN2A.p61
CDKN2A.p69
CDKN2A.p69
ERBB2.p755
ERBB2.p755
```
outputs
```
#!text
CDKN2A.p61      CDKN2A.p.75    9       21971175-21971176-21971177      21971176-21971177-21971178      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p61      CDKN2A.p.76    9       21971175-21971176-21971177      21971173-21971174-21971175      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      CDKN2A.p.55    9       21971194-21971195-21971196      21971193-21971194-21971195      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS],CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      CDKN2A.p.84    9       21971151-21971152-21971153      21971149-21971150-21971151      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      CDKN2A.p.83    9       21971151-21971152-21971153      21971152-21971153-21971154      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      CDKN2A.p.54    9       21971194-21971195-21971196      21971196-21971197-21971198      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      CDKN2A.p.55    9       21971194-21971195-21971196      21971193-21971194-21971195      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS],CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69      CDKN2A.p.84    9       21971151-21971152-21971153      21971149-21971150-21971151      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      CDKN2A.p.83    9       21971151-21971152-21971153      21971152-21971153-21971154      CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69      CDKN2A.p.54    9       21971194-21971195-21971196      21971196-21971197-21971198      CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
ERBB2.p755      ERBB2.p.725   17      37880219-37880220-37880221      37880219-37880220-37880221      CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
ERBB2.p755      ERBB2.p.785   17      37881024-37881025-37881026      37881024-37881025-37881026      CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
ERBB2.p755      ERBB2.p.725   17      37880219-37880220-37880221      37880219-37880220-37880221      CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS],CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
ERBB2.p755      ERBB2.p.785   17      37881024-37881025-37881026      37881024-37881025-37881026      CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS],CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
```
The third column indicates the potential transcript usage for the alternative identifier. Each transcript usage is denoted by <listing transcript>/<actual transcript>. Different potential choices are separated by ','.

#### annotate mutations from genomic locations

This is the forward annotation

```
#!bash
revan anno --ccds CCDS.current.txt --ref hs37d5.fa -i 'chr3:178936091.G>A'
```
outputs
```
#!text
chr3:178936091.G>A      3       178936091-178936092-178936093   CCDS43171.1
     PIK3CA (+, coding)      3:G178936091A/c.1633>/p.E545K   .
```

### Technical notes

RevAn follows in full the HGVS nomenclature while annotating protein level mutation identifiers. For example, a out-of-phase, in frame insertion, `ACIN1:c.1930_1931insATTCAC` will be annotated with `p.S643_R644insHS` rather than `R644delinsHSR`. Protein level mutation will be generated as if no nucleotide mutation information exists.

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.