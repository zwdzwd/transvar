**TransVar** is a reverse annotator for inferring genomic characterization(s) of mutations (e.g., ```chr3:178936091 G=>A```) from transcript-dependent annotation(s) (e.g., ```PIK3CA p.E545K``` or ```PIK3CA c.1633G>A```). It is designed for resolving ambiguous mutation annotations arising from differential transcript usage. TransVar has the following features:

 + supports HGVS nomenclature
 + supports both left-alignment and right-alignment convention in reporting indels.
 + supports annotation of a region based on a transcript dependent characterization
 + supports single nucleotide variation (SNV), insertions and deletions (indels) and block substitutions
 + supports mutations at both coding region and intronic/UTR regions
 + supports transcript annotation from commonly-used databases such as Ensembl, NCBI RefSeq and GENCODE etc
 + supports UniProt protein id as transcript id
 + supports GRCh36, 37, 38
 + functionality of forward annotation.

--------

[TOC]

--------

### Quick start

download program, unzip and cd to directory
```
#!bash
./transvar config --download_hg19_anno
./transvar revanno --custom -i 'PIK3CA.p.E545K'
```

Note that to use TransVar in full, one need to link to a (samtools faidx indexed) reference assembly. One can set the default location (e.g., ./hg19.fa) for the reference assembly via
```
#!bash
transvar config -k reference -v ./hg19.fa -s hg19
```
Please see Install section for detailed instruction.

### Download and install

#### dependency

Basic functionalities requires just Python 2.7. Some additional annotation also depends on the pysam library.

#### program

current stable version: [version 1.28](https://bitbucket.org/wanding/transvar/get/v1.28.zip)

#### reference genome assembly
For most annotation tasks, TransVar requires a samtools faidx indexed reference genome in fasta format, which is available at, e.g., [UCSC ftp](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/).
Once downloaded and indexed, the reference can be used through the "--reference" option followed by the fasta filename.

To set the default location of reference to ./hg19.fa,
```
#!bash
transvar config -k reference -v ./hg19.fa -s hg19
```
will create in transvar.cfg an entry
```
#!text
[hg19]
reference = hg19.fa
```
so that there is no need to specify the location of reference on subsequent usages.

#### transcript annotations

TransVar provides automatic download of transcript annotations specific to common versions of human genome.
```
#!bash
transvar config --download_hg19_anno
```
will automatically download annotation from Ensembl, RefSeq etc. to [install dir]/transvar.download directory or your local ~/.transvar.download if the installation directory is inaccessible.
See "transvar config -h" for downloading more versions.
These will also create default mappings under the corresponding reference version section of transvar.cfg like
```
#!text
[hg19]
ucsc = /home/wzhou1/download/hg19.ucsc.txt.gz
```

##### Annotating nonprotein coding genomic features
In annotating non-protein-coding genomic sequences such as lincRNA, pseudogenes etc. TransVar requires a tab-indexed transcript database.
For example,
```
#!bash
zgrep -v '^#' download/hg19.gencode.gtf.gz | sort -k1,1 -k4,4n | bgzip > download/hg19.gencode.sorted.gtf.gz
tabix -p gff download/hg19.gencode.sorted.gtf.gz
```

### Usage

#### specify transcript annotation

The following table summarize the transcript annotations supported by TransVar as well as the option(s) to use each database in the annotation. TransVar can take any one or any combination(s) of these transcript annotations as long as these annotations are based on the same version of reference assembly.

 | Database | Format | Default  | Non-default     |
 |:-------|:--------|:---------|:----|
 | CCDS | CCDS flat text | `--ccds` | `--ccds CCDS.current.txt` |
 | UCSC | UCSC RefGene | `--ucsc` | `--ucsc2 hg19.ucsc.txt` |
 | Ensembl | Ensembl GTF | `--ensembl`  | `--ensembl GRCh37.75.gtf.gz`  |
 | RefSeq | RefSeq GFF3 | `--refseq`  | `--refseq GRCh37.p13.gff3.gz`   |
 | AceView | AceView GFF | `--aceview` | `--aceview AceView.ncbi37.gff.gz`  |
 | GENCODE | GENCODE GTF | `--gencode` | `--gencode gencode.v19.gtf.gz`  |
 | knownGene | knownGene table | `-kg` | `--kg kg.gz --alias kgAlias.gz` |
 | custom | custom table | `--custom` | `--custom hg19.map` |

If one download transcripts through "transvar config", TransVar would use the downloaded definition automatically (by setting the default configuration file). For example, "--ccds" would look for the downloaded CCDS definition. One can specify non-default annotation by appending a path to the option ("--ccds CCDS.current.txt").
To set the default annotation of a particular reference version,
```
#!bash
transvar config -k ccds -v CCDS.current.txt -s hg19
```
The configuration file is located either at the "[install dir]/transvar.cfg" or "~/.transvar.cfg" if the installation directory is inaccessible.

---

#### specify reference assembly

TransVar provide native support for switching between reference assemblies. Each reference assembly is represented in a section such as
```
#!text
[DEFAULT]
refversion = hg19

[hg19]
refseq = /home/wzhou1/transvar/download/hg19.refseq.gff.gz
ccds = /home/wzhou1/transvar/download/hg19.ccds.txt
ensembl = /home/wzhou1/transvar/download/hg19.ensembl.gtf.gz
reference = /projects/database/reference/hg19.fa

[hg38]
refseq = /home/wzhou1/transvar/download/hg38.refseq.gff.gz
gencode = /home/wzhou1/transvar/download/hg38.gencode.gtf.gz
ucsc = /home/wzhou1/transvar/download/hg38.ucsc.txt.gz
reference = /home/wzhou1/reference/hg38.fa
```
The "refversion" key specify the default reference version ("hg19" in the above example).

To add a new version and specify the location of some transcript annotation
```
#!bash
transvar config -k ccds -v ccds.myhg.txt -s myhg
```
Will create in transvar.cfg a section like
```
#!text
[myhg]
ccds = ccds.myhg.txt
```

To switch to a version on the fly, one could use the "--refversion" option, e.g.,
```
#!bash
transvar revanno -i 'PIK3CA:p.E545K' --ucsc --refversion hg38
```

To change the default reference,
```
#!bash
transvar config -k refversion -v hg38
```

#### batch processing

For all mutation types, one can batch process a list of mutation identifiers with optional transcript id to constraint the search. Take SNV for example,
```
#!bash
transvar revanno -l example/input_table -g 1 -m 4 -t 2 --ensembl -o 2,3,5
```
As suggested by the command, TransVar takes as input the 1st column as gene and 4th column as identifier. The 2nd column will be used as the transcript id from Ensembl to constrain the alternative identifier search. The 2nd, 3rd and 5th columns are chosen to be output as a validation of TransVar's performance.

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


#### reverse annotation of protein sites

To use uniprot id as protein name, one must first download the uniprot id map by
```
#!bash
transvar config --download_idmap
```

Then one could use protein id instead of gene name by applying the `--uniprot` option to TransVar. For example,

```
#!bash
transvar revanno --ccds -i 'Q5VUM1:47' --uniprot
```
```
#!text
Q5VUM1:47   6   71289191-71289193   CCDS4972.1  C6ORF57 (+, coding)
    6:g.71289191_71289193/c.139_141/p.47    PRefSeq=S;NRefSeq=TCC;RefSeq=TCC
```
TransVar use a keyword extension `ref` in `Q5VUM1:p.47refS` to differentiate from the synonymous mutation `Q5VUM1:p.47S`. The former notation specifies that the reference protein sequence is `S` while the later specifies the target protein sequence is `S`.

---

#### reverse annotation of protein motif

For example, one can find the genomic location of a DRY motif in protein P28222 by issuing the following command,
```
#!bash
transvar revanno -i 'P28222.p.146_148refDRY' --uniprot --ccds
```
```
#!text
P28222.p.146_148refDRY   6   78172677-78172685    CCDS4986.1
    HTR1B (-, coding)  6:g.78172677_78172685/c.436_444/p.146_148
    PRefSeq=DRY;NRefSeq=GACCGCTAC;RefSeq=GTAGCGGTC
```
One can also use wildcard `x` (lowercase) in the motif.
```
#!bash
transvar revanno -i 'HTR1B.p.365_369refNPxxY' --ccds
```
```
#!text
HTR1B.p.365_369refNPxxY 6   78172014-78172028    CCDS4986.1
    HTR1B (-, coding)   6:g.78172014_78172028/c.1093_1107/p.365_369
    PRefSeq=NPIIY;NRefSeq=AAC..TAT;RefSeq=ATA..GTT
```

---
#### reverse annotation of single amino acid substitution
Mutation formats acceptable in TransVar are ```PIK3CA:E545K``` or ```PIK3CA:p.E545K```, or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. TransVar takes native HGVS format inputs and outputs. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.

```
#!bash
transvar revanno -i PIK3CA:E545K --ensembl
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
#!bash
transvar revanno -i ACSL4:p.R133R --ccds
```
```
#!text
ACSL4   23:108926078    108926078       c.399C>T        p.R133R Missense
    X       108926078-108926079-108926080   CCDS14548.1
    ACSL4 (-, coding)       X:G108926078A/c.399C>T/p.R133R
    CddMuts=X:G108926078T,X:G108926078C,X:G108926078A;NCodonSeq=CGC;NCddSeqs=AGG,AGA,CGA,CGC,CGG,CGT
```
In those cases, TransVar prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

If one wishes, he/she can also turn on the `--dbsnp` option to show further potential match to the dbSNP database. (This requires 1) pysam library installed in python path; 2) dbsnp library downloaded via `transvar config --download_hg19_dbsnp`)
```
#!bash
$ transvar revanno -i 'A1CF:p.A309A' --ccds --dbsnp
```
```
#!text
A1CF:p.A309A    10    52576004-52576005-52576006    CCDS7243.1
    A1CF (-, coding, exon 7)    10:g.52576004T>G/c.927A>C/p.A309A
    NCodonSeq=GCA;NCddSeqs=GCC,GCG,GCT;CddSNVMuts=10:g.52576004T>C,10:g.52576004T>A;
	DBSNP=rs201831949(10:52576004T>G)
```
Note that in order to use dbSNP, one must download the dbSNP database through `transvar config --download_hg19_dbsnp`, or by configure the `dbsnp` slot in the configure file via `transvar config -k dbsnp -v [path to dbSNP VCF]`. dbSNP file must be tabix indexed.


---
#### reverse annotation of single nucleotide variation (SNV)

TransVar infers nucleotide mutation through ```PIK3CA:1633G>A``` or ```PIK3CA:c.1633G>A```. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.
```
#!bash
transvar revanno --ccds -i 'PIK3CA:c.1633G>A'
```
outputs
```
#!text
PIK3CA:c.1633G>A     3     178936091-178936092-178936093   CCDS43171.1
	PIK3CA (+ coding)  3:G178936091A/c.1633G>A/p.E545K NCodonSeq=GAG;NAltCodonSeq=AAG
```
---

#### reverse annotation of nucleotide insertion
An insertion may result in: 1) a pure insertion of amino acids; 2) a block substitution of amino acids, when insertion occur after 1st or 2nd base in a codon; or 3) a frame-shift. Following HGVS nomenclature, TransVar labels the first different amino acid and the length of the peptide util stop codon, assuming no change in the splicing.

Example: to annotate an **in-frame, in-phase insertion**,
```
#!bash
transvar revanno --ccds -i 'ACIN1:c.1932_1933insATTCAC'
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
transvar revanno --ccds -i 'ACIN1:c.1930_1931insATTCAC'
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
transvar revanno --ccds -i 'AAAS:c.1225_1226insG'
```
```
#!text
AAAS:c.1225_1226insG    12   12:53702089-53702090    CCDS8856.1    AAAS (-, coding)
    12:53702089_53702090insC/c.1225_1226insG/p.E409Gfs*17   NatInsSeq=G;RefInsSeq=C
```

Example: to annotate an **intronic insertion**,
```
#!bash
transvar revanno --ccds -i 'ADAM33:c.991-3_991-2insC'
```
```
#!text
ADAM33:c.991-3_991-2insC   20   20:3654141-3654142-3654143-(3654145)-(ins)-(3654146)
    CCDS13058.1     ADAM33 (-, intronic)     20:3654145_3654146insG/c.991-3_991-2insC/.
    RefInsSeq=G;NatInsSeq=C
```
In the case of intronic insertions, amino acid identifier is not applicable, represented in a `.`.

---

#### reverse annotation of nucleotide deletion
Similar to insertions, deletion can be in-frame or frame-shift. The consequence of deletion to amino acid sequence may appear a simple deletion or a block substitution (in the case where in-frame deletion is out of phase, i.e., partially delete codons).

Example: to annotate an **in-frame deletion**,
```
#!bash
transvar revanno --ccds -i 'A4GNT:c.694_696delTTG'
```
```
#!text
A4GNT:c.694_696delTTG   3   137843433-137843435     CCDS3097.1    A4GNT (- coding)
    3:137843433_137843435del/c.694_696del/p.L232del RefDelSeq=CAA;NatDelSeq=TTG
```

Example: to annotate a **in-frame, out-of-phase deletion**,
```
#!bash
transvar revanno --ccds -i 'ABHD15.c.431_433delGTG'
```
```
#!text
ABHD15.c.431_433delGTG  17   27893552-27893554    CCDS32602.1     ABHD15 (- coding)
    17:27893552_27893554del/c.431_433del/p.C144_V145delinsF RefDelSeq=CAC;NatDelSeq=GTG
```

Example: to annotate a **frame-shift deletion**,
```
#!bash
transvar revanno --ccds -i 'AADACL3.c.374delG'
```
```
#!text
AADACL3.c.374delG    1    12785494-12785494    CCDS41252.1   AADACL3 (+ coding)
    1:12785494_12785494del/c.374del/p.C125Ffs*17    RefDelSeq=G;NatDelSeq=G
```

Example: to annotate a **deletion that span from intronic to coding region**,
```
#!bash
transvar revanno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'
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
transvar revanno --ccds -i 'A1CF:c.508_509CC>TT'
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
transvar revanno --ccds -i 'CSRNP1.c.1212_1224>GGAGGAGGAA'
```
```
#!text
CSRNP1.c.1212_1224>GGAGGAGGAA   3    39185092-39185104   CCDS2682.1   CSRNP1 (-, coding)
    3:39185092_39185104TTCCTCCTCCTCC>TTCCTCCTCC/c.1212_1224GGAGGAGGAGGAA>GGAGGAGGAA/p.E408del .
```

Example: to annotate a block substitution in **intronic region**,
```
#!bash
transvar revanno --ccds -i 'A1CF:c.1460+2_1460+3TG>CC'
```
```
#!text
A1CF:c.1460+2_1460+3TG>CC    10    52570797-52570798   CCDS7241.1    A1CF (-, intronic)
    10:52570797_52570798CA>GG/c.1460+2_1460+3TG>CC/.        .
```

When block substitution occurs **across splice site**, TransVar put a tag in the info fields and does not predict amino acid change.
```
#!bash
transvar revanno --ccds -i 'A1CF:c.1459_1460+3ATGTG>CC'
```
```
#!text
A1CF:c.1459_1460+3ATGTG>CC    10   52570797-52570801   CCDS7241.1   A1CF (-, coding;intronic)
	10:52570797_52570801CACAT>GG/c.1459_1460+3ATGTG>CC/.    CrossSplitSite
```

---

#### reverse annotation of nucleotide duplication

Duplication can be thought of as insertion where the inserted sequence is identical to the sequence flanking the breakpoint.
Similar to insertion, the annotation of duplication assumes no change in splicing.

Example: to annotate a duplication coding region,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1669_1674dup'
```
```
#!text
CHD7:c.1669_1674dup    8    61693562-61693567 (dup) CCDS47865.1     CHD7 (+, Coding)
    8:61693562-61693567dupTCCCCG/c.1669_1674dup/p.S557_P558dupSP
    RefDupSeq=TCCCCG;NatDupSeq=TCCCCG
```

Example: a duplication on the nucleotide level may lead to frame-shift or block substitution on the amino acid level,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1668_1669dup'
```
```
#!text
CHD7:c.1668_1669dup    8    61693561-61693562 (dup) CCDS47865.1     CHD7 (+, Coding)
    8:61693561-61693562dupTT/c.1668_1669dup/p.S557Ffs*8   RefDupSeq=TT;NatDupSeq=TT
```

Example: to annotate a duplication in intronic region,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1666-5_1666-3dup'
```
```
#!text
CHD7:c.1666-5_1666-3dup 8   61693554-61693556 (dup) CCDS47865.1   CHD7 (+, Intronic)
    8:61693554-61693556dupCTC/c.1666-5_1666-3dup/.  RefDupSeq=CTC;NatDupSeq=CTC
```


---

#### reverse annotation of amino acid insertion

```
#!bash
transvar revanno --ccds -i 'AATK.p.P1331_A1332insTP'
```
```
#!text
AATK    c.3993_3994insACGCCC   p.P1331_A1332insTP   17:79093270-79093271
    17    79093268-79093273 (insertion)   CCDS45807.1     AATK (-, coding)
    17:79093268-79093273ins6/c.3991-3996ins6/p.P1331_A1332insTP     Uncertain
```

#### reverse annotation of amino acid deletion
```
#!bash
transvar revanno --ccds -i 'AADACL4.p.W263_I267delWRDAI'
```
```
#!text
AADACL4   c.788_802del15  p.W263_I267delWRDAI   1:12726310-12726324
     1       12726309-12726323 (deletion)    CCDS30590.1     AADACL4 (+, coding)
     1:12726309-12726323/c.787-801/p.W263_I267delWRDAI     Uncertain
```

#### reverse annotation of amino acid block substitution
```
#!bash
transvar revanno --ccds -i 'ABCC3:p.Y556_V557delinsRRR'
```
```
#!text
ABCC3:p.Y556_V557delinsRRR   17   48745254-48745259 (block substitution)  CCDS32681.1
    ABCC3 (+, coding)
    17:48745254-48745259TACGTG>AGGAGGAGG/c.1666-1671TACGTG>AGGAGGAGG/p.Y556_V557delinsRRR
    CddNatAlt=AGG/AGA/CGA/CGC/CGG/CGT+AGG/AGA/CGA/CGC/CGG/CGT+AGG/AGA/CGA/CGC/CGG/CGT;Uncertain
```

#### reverse annotation of amino acid frame-shift

```
#!bash
transvar revanno --ccds -i 'A1BG.p.G132fs*2'
```
```
#!text
A1BG.p.G132fs*2 19      58863866-58863867-58863868      CCDS12976.1     A1BG (-, coding)
    19:58863860-58863868/c.394-402/p.G132fs*2       RoughEstimateFromFrameShift
```

---

#### search alternative codon identifiers

An identifier is regarded as an alternative if the underlying codon overlap with the one from the original identifier.
Example: to search alternative identifiers of CDKN2A.p.58 (without knowing reference allele),
```
#!bash
transvar codonsearch --ccds -i CDKN2A.p.58
```
```
#!text
CDKN2A.p.58    CDKN2A.p.72   9    21971184-21971185-21971186   21971185-21971186-21971187
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p.58    CDKN2A.p.73   9    21971184-21971185-21971186   21971182-21971183-21971184
    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
```
The pair of transcript id listed corresponds to the transcripts based on which, the original and alternative identifiers are defined. Multiple pairs of transcript definitions are appended following a `,`.

Example: to search alternative identifiers of DHODH:G152R (knowing reference allele `G`, alternative allele here will be ignored),
```
#!bash
transvar codonsearch -i DHODH:G152R --refseq
```
outputs
```
#!text
DHODH:G152R     DHODH.p.G9    16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq],NM_001361.4[RefSeq]/XM_005255829.1[RefSeq]
DHODH:G152R     DHODH.p.G124  16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq],NM_001361.4[RefSeq]/XM_005255827.1[RefSeq]
DHODH:G152R     DHODH.p.G16   16      72050942-72050943-72050944      72050942-72050943-72050944      NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq],NM_001361.4[RefSeq]/XM_005255828.1[RefSeq]
```
TransVar outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

Example: to run `transvar codonsearch` to **batch process** a list of mutation identifiers.
```
#!bash
transvar codonsearch -l example/input.txt --ccds -m 1 -o 1
```
Example input table
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

---
#### infer potential codon identity

Example: to check if MET.p1010 and MET.p992 may be refering to one mutation due to different usage of transcripts,
```
#!bash
transvar codonsearch --refseq -i MET.p1010
```
gives
```
#!text
MET.p.1010      MET.p.1047      7       116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/XM_005250353.1[RefSeq]
MET.p.1010      MET.p.580       7       116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/XM_005250354.1[RefSeq]
MET.p.1010      MET.p.991       7       116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/NM_001127500.1[RefSeq]
MET.p.1010      MET.p.1029      7       116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/XM_005250353.1[RefSeq]
MET.p.1010      MET.p.543       7       116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/XM_005250354.1[RefSeq]
MET.p.1010      MET.p.992       7       116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/NM_000245.2[RefSeq]
MET.p.1010      MET.p.973       7       116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/NM_000245.2[RefSeq]
MET.p.1010      MET.p.1028      7       116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/NM_001127500.1[RefSeq]
MET.p.1010      MET.p.562       7       116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/XM_005250354.1[RefSeq]
```
Since MET.p.992 is in the list, the two identifiers might be due to the same genomic mutation.

#### annotate SNP from genomic locations

This is the forward annotation

```
#!bash
transvar anno --ccds -i 'chr3:178936091.G>A'
```
outputs
```
#!text
chr3:178936091.G>A   3   178936091-178936092-178936093   CCDS43171.1
     PIK3CA (+, coding)      3:G178936091A/c.1633>/p.E545K   .
```

Another example:
```
#!bash
transvar anno -i "9:135782704C>G" --ccds
```
outputs
```
#!text
9:135782704C>G  9    135782704    CCDS6956.1    TSC1 (-, coding)
    9:g.135782704C>G/c.1317G>C/p.L439L
    CodonPos=135782704-135782705-135782706;NCodonSeq=CTG
9:135782704C>G  9    135782704    CCDS55350.1   TSC1 (-, coding)
    9:g.135782704C>G/c.1164G>C/p.L388L
    CodonPos=135782704-135782705-135782706;NCodonSeq=CTG
```


#### annotate a short genomic region

To annotate a short genomic region in a gene,
```
#!bash
transvar anno --ccds -i 'chr3:g.178936091_178936192'
```
outputs
```
#!text
chr3:g.178936091_178936192   3    178936091-178936192    CCDS43171.1
    PIK3CA (+, coding,intronic)   3:g.178936091_178936192/c.1633_1664+70/p.E545_R555
    BEGCodon=178936091-178936092-178936093;ENDCodon=178936121-178936122-178936984
```
	
Results indicates the beginning position is at coding region while ending position is at intronic region (c.1633_1664+70).

For intergenic sites, TransVar also reports the identity and distance to the gene upstream and downstream. For example, `chr6:116991832` is simply annotated as intergenic in the original annotation. TransVar reveals that it is 1,875 bp downstream to ZUFSP and 10,518 bp upstream to KPNA5 showing a vicinity to the gene ZUFSP. There is no limit in the reported distance. If a site is at the end of the chromosome, TransVar is able to report the distance to the telomere.

#### annotate a long genomic region
[back to top](#top)
```
#!bash
transvar anno -i '9:g.133750356_137990357' --ccds
```
outputs
```
#!text
9:g.133750356_137990357 9     133750356-137990357     BEG=CCDS35165.1,END=CCDS6986.1
    4,240,002 bp covering 53 genes  9:g.133750356_137990357/./.
    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1244A>/p.H415;
	ENDreg=OLFM1 (+, intronic);ENDid=9:g.137990357C>/c.622+6C>/.
9:g.133750356_137990357 9     133750356-137990357     BEG=CCDS35166.1,END=CCDS6986.1
    4,240,002 bp covering 53 genes  9:g.133750356_137990357/./.
    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1187A>/p.H396;
	ENDreg=OLFM1 (+, intronic);ENDid=9:g.137990357C>/c.622+6C>/.
```
The result indicates that the region span 53 genes. The beginning of the region resides in the coding sequence of ABL1, c.1187A and the ending region resides in the intronic region of OLFM1, c.622+6C. 2 different usage of transcripts in annotating the starting position is represented in two lines, each line corresponding to a combination of transcript usage.
This annotation not only shows the coverage of the region, also reveals the fine structure of the boundary.

In another example, where the ending position exceeds the length of the chromosome, TransVar truncates the region and outputs upstream and downstream information of the ending position.
```
#!bash
transvar anno -i '9:g.133750356_1337503570' --ccds
```
outputs
```
#!text
9:g.133750356_1337503570    9    133750356-141213431    BEG=CCDS35165.1,END=.
    7,463,076 bp covering 137 genes 9:g.133750356_141213431/./.
    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1244A>/p.H415;
	ENDreg=Noncoding (up: 484,026 bp to EHMT1, down: 0 bp to 3-telomere);ENDid=././.
9:g.133750356_1337503570    9    133750356-141213431    BEG=CCDS35166.1,END=.
    7,463,076 bp covering 137 genes 9:g.133750356_141213431/./.
    BEGreg=ABL1 (+, coding);BEGid=9:g.133750356A>/c.1187A>/p.H396;
	ENDreg=Noncoding (up: 484,026 bp to EHMT1, down: 0 bp to 3-telomere);ENDid=././.
```

#### annotate a deletion from genomic location
[back to top](#top)

A frameshift deletion
```
#!bash
transvar anno -i "chr2:234183368_234183380del" --ccds
```
outputs
```
#!text
chr2:234183368_234183380del     chr2    234183368-234183380     CCDS2502.2
    ATG16L1 (+, Coding)
    chr2:g.234183368_234183380del13/c.841_853del13/p.T281Lfs*5
    LEFTALNG=g.234183367_234183379del13;
	UALNG=g.234183368_234183380del13;
	LEFTALNC=c.840_852del13;
	UALNC=c.841_853del13;
	REG=Exonic_8
```
Note the difference between left-aligned identifier and the right aligned identifier.

An in-frame deletion
```
#!bash
transvar anno -i "chr2:234183368_234183379del" --ccds
```
outputs
```
#!text
chr2:234183368_234183379del     chr2    234183368-234183379     CCDS2502.2
    ATG16L1 (+, Coding)     chr2:g.234183368_234183379del12/c.841_852del12/p.T281_G284del
	BEGCodon=234183368/234183369/234183370;
	ENDCodon=234183377/234183378/234183379;REG=Exonic_8
```

Another example
```
#!bash
transvar anno --ccds -i '12:53703425_53703427del'
```
outputs
```
#!text
12:53703425_53703427del chr12   53703425-53703427       CCDS8856.1
    AAAS (-, Coding)
    chr12:g.53703427_53703429delCCC/c.769_771delGGG/p.257delG
    LEFTALNG=g.53703424_53703426delCCC;
	UALNG=g.53703425_53703427delCCC;
	LEFTALNC=c.766_768delGGG;
	UALNC=c.768_770delGGG;
	REG=Exonic_8;
	LEFTALNP=p.256delG;
	UALNP=p.256delG
```
Note the difference between left and right-aligned identifiers on both protein level and cDNA level.

An in-frame out-of-phase deletion
```
#!bash
transvar anno -i "chr2:234183372_234183383del" --ccds
```
outputs
```
#!text
chr2:234183372_234183383del     chr2    234183372-234183383     CCDS2502.2
    ATG16L1 (+, Coding)     chr2:g.234183372_234183383del12/c.845_856del12/p.H282_G286delinsR
    BEGCodon=234183371/234183372/234183373;
	ENDCodon=234183383/234183384/234183385;REG=Exonic_8
```

#### annotate an insertion from genomic location

An in-frame insertion of three nucleotides
```
#!bash
transvar anno -i '2:69741762insTGC' --ccds
```
outputs
```
#!text
2:69741762insTGC        chr2    69741762        CCDS1893.2
    AAK1 (-, Exonic_12)
    chr2:g.69741762insTGC/c.1616_1617insGCA/p.Q546_L547insQ
	LEFTALNP=p.Y532_Q533insQ
```
Note the proper right-alignment of protein level insertion Q. The left-aligned identifier is also given in the `LEFTALN` field.

A frame-shift insertion of two nucleotides
```
#!bash
transvar anno -i '7:121753754insCA' --ccds
```
outputs
```
#!text
7:121753754insCA    chr7    121753754      CCDS5783.1
    AASS (-, Exonic_9)
    chr7:g.121753754insCA/c.1063_1064insTG/p.I355Mfs*10  .
```

```
#!bash
transvar anno -i '17:79093270insGGGCGT' --ccds
```
outputs
```
17:79093270insGGGCGT    chr17   79093270        CCDS45807.1
    AATK (-, Exonic_13)
    chr17:g.79093273insCGTGGG/c.3993_3994insACGCCC/p.P1331_A1332insTP
    LEFTALNG=g.79093270insGGGCGT;UALNG=g.79093270insGGGCGT;
	LEFTALNC=c.3976_3977insCGCCCA;UALNC=c.3993_3994insACGCCC;
	LEFTALNP=p.A1326_P1327insPT;UALNP=p.P1331_A1332insTP
```
Notice the difference in the inserted sequence when left-alignment and right-alignment conventions are followed.

A frame-shift insertion of one nucleotides in a homopolymer
```
#!bash
transvar anno -i '7:117230474insA' --ccds
```
outputs
```
#!text
7:117230474insA chr7    117230474       CCDS5773.1    CFTR (+, Exonic_13)
    chr7:g.117230474insA/c.1752_1753insA/p.E585Rfs*4
    LEFTALNC=c.1747_1748insA
```
Notice the right alignment of cDNA level insertion and the left alignment reported as additional information.

A in-frame, in-phase insertion
```
#!bash
transvar anno -i '12:109702119insACC' --ccds
```
```
#!text
12:109702119insACC   chr12   109702119    CCDS31898.1
    ACACB (+, Exonic_49)
    chr12:g.109702119insACC/c.6870_6871insACC/p.Y2290_H2291insT   .
```

#### annotate block substitution from genomic locations

A block-substitution that results in a frameshift.
```
#!bash
transvar anno -i 'chr10:g.27329002_27329002A>AT' --ccds
```
```
#!text
chr10:g.27329002_27329002A>AT   chr10   27329002-27329002
    CCDS41499.1     ANKRD26 (-, CDS_21)
    chr10:g.27329002_27329002A>AT/c.2267_2267T>AT/p.M756Nfs*6
    BEGCodon=27329001-27329002-27329003;ENDCodon=27329001-27329002-27329003
```

A block-substitution that is in-frame,
```
#!bash
transvar anno -i 'chr10:g.52595929_52595930GG>AA' --ccds
```
```
#!text
chr10:g.52595929_52595930GG>AA  chr10   52595929-52595930
    CCDS7241.1      A1CF (-, CDS_4) chr10:g.52595929_52595930GG>AA/c.508_509CC>TT/p.P170L
    BEGCodon=52595928-52595929-52595930;ENDCodon=52595928-52595929-52595930
```

### FAQ

+ I got 'GeneNotRecognized', what's wrong?

Most likely you forgot to specify a transcipt definition such as `--ccds` or `--ensembl`. Sometimes there are non-canonical names for genes, this can be fixed through the `--alias` option and specify an alias table. TransVar comes with alias table from UCSC knownGene.

### Technical notes

TransVar follows in full the HGVS nomenclature while annotating protein level mutation identifiers. For example, a out-of-phase, in frame insertion, `ACIN1:c.1930_1931insATTCAC` will be annotated with `p.S643_R644insHS` rather than `R644delinsHSR`. Protein level mutation will be generated as if no nucleotide mutation information exists.

## Future work

 + add cytoband annotation
 + forward annotation of splice site
 + forward annotation of non-coding RNA from GENCODE
 + forward annotation of binding sites
 + forward annotation of structural variation breakpoints
 + support uncertain insertion such as '4424_4425ins80'
 + bundle Pysam

## Bug report and feature request

If you find any bug (very likely due to the complexity of genomics:-)) or you wish any feature, please direct to Wanding Zhou <wzhou1@mdanderson.org>. Thank you.

## Reference

We are working on an application note on this topic :-).

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.
