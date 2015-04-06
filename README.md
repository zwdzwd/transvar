**TransVar** is a reverse annotator for inferring genomic characterization(s) of mutations (e.g., ```chr3:178936091 G=>A```) from transcript-dependent annotation(s) (e.g., ```PIK3CA p.E545K``` or ```PIK3CA c.1633G>A```). It is designed for resolving ambiguous mutation annotations arising from differential transcript usage. TransVar keeps awareness of the underlying unknown transcript structure (exon boundary, reference amino acid/base) while performing reverse annotation.
TransVar has the following features:

 + supports HGVS nomenclature
 + supports both left-alignment and right-alignment convention in reporting indels and duplications.
 + supports annotation of a region based on a transcript-dependent characterization
 + supports noncoding RNA annotation
 + supports single nucleotide variation (SNV), insertions and deletions (indels) and block substitutions
 + supports mutations at both coding region and intronic/UTR regions
 + supports transcript annotation from commonly-used databases such as Ensembl, NCBI RefSeq and GENCODE etc
 + supports UniProt protein id as transcript id
 + supports GRCh36, 37, 38 (human),  GRCm38 (mouse), NCBIM37 (mouse)
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

current stable version: [version 1.36](https://bitbucket.org/wanding/transvar/get/v1.36.zip)

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

<!-- ##### Annotating nonprotein coding genomic features -->
<!-- In annotating non-protein-coding genomic sequences such as lncRNA, pseudogenes etc. TransVar requires a tab-indexed transcript database. -->
<!-- For example, -->
<!-- ``` -->
<!-- #!bash -->
<!-- zgrep -v '^#' download/hg19.gencode.gtf.gz | sort -k1,1 -k4,4n | bgzip > download/hg19.gencode.sorted.gtf.gz -->
<!-- tabix -p gff download/hg19.gencode.sorted.gtf.gz -->
<!-- ``` -->

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

#### view current configuration

One can read the transvar.cfg file for the information. Alternatively one may run
```
#!bash
transvar current
```
which returns information about the setup regarding to the current reference selection, including the location of the reference file and database file.
```
#!text
Current reference version: mm10
reference: /home/wzhou/genomes_link/mm10/mm10.fa
Available databases:
refseq: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.refseq.gff.gz
ccds: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.ccds.txt
ensembl: /home/wzhou/tools/transvar/transvar/transvar.download/mm10.ensembl.gtf.gz
```

#### batch processing

For all mutation types, one can batch process a list of mutation identifiers with optional transcript id to constraint the search. Take SNV for example,
```
#!bash
transvar revanno -l example/input_table -g 1 -m 5 -t 2 --ensembl -o 2,3,4
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
input   chrm    pos     transcript      gene    strand  coordinates(gDNA/cDNA/protein)  region  info
ENST00000286744|15:84442328|c.243G>A    chr15   84442326-84442327-84442328      ENST00000286744 ADAMTSL3        +       chr15:g.84442327G>A/c.242G>A/p.W81*     cds_in_exon_4   reference_codon=TGG;candidate_codons=TAA,TAG,TGA;candidate_snv_variants=chr15:g.84442328G>A;candidate_mnv_variants=chr15:g.84442327_84442328GG>AA
ENST00000286744|15:84442326|c.241T>C    chr15   84442326-84442327-84442328      ENST00000286744 ADAMTSL3        +       chr15:g.84442326T>A/c.241T>A/p.W81R     cds_in_exon_4   reference_codon=TGG;candidate_codons=AGG,AGA,CGA,CGC,CGG,CGT;candidate_snv_variants=chr15:g.84442326T>C;candidate_mnv_variants=chr15:g.84442326_84442328TGG>AGA,chr15:g.84442326_84442328TGG>CGA,chr15:g.84442326_84442328TGG>CGC,chr15:g.84442326_84442328TGG>CGT
ENST00000369038|1:150530513|c.2270G>A   chr1    150530512-150530513-150530514   ENST00000369038 ADAMTSL4        +       chr1:g.150530513G>A/c.2270G>A/p.G757D   cds_in_exon_12  reference_codon=GGT;candidate_codons=GAC,GAT;candidate_mnv_variants=chr1:g.150530513_150530514GT>AC
ENST00000338316|5:7802364|c.2662G>A     chr5    7802364-7802365-7802366 ENST00000338316 ADCY2   +       chr5:g.7802364G>A/c.2662G>A/p.V888I     cds_in_exon_21  reference_codon=GTC;candidate_codons=ATC,ATA,ATT;candidate_mnv_variants=chr5:g.7802364_7802366GTC>ATA,chr5:g.7802364_7802366GTC>ATT
ENST00000338316|5:7802365|c.2663T>C     chr5    7802364-7802365-7802366 ENST00000338316 ADCY2   +       chr5:g.7802365T>C/c.2663T>C/p.V888A     cds_in_exon_21  reference_codon=GTC;candidate_codons=GCA,GCC,GCG,GCT;candidate_mnv_variants=chr5:g.7802365_7802366TC>CA,chr5:g.7802365_7802366TC>CG,chr5:g.7802365_7802366TC>CT
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
Q5VUM1:47       chr6    71289191-71289193       CCDS4972.1      C6ORF57 +
  chr6:g.71289191_71289193/c.139_141/p.47S        cds_in_exon_2
  protein_sequence=S;cDNA_sequence=TCC;gDNA_sequence=TCC
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
P28222.p.146_148refDRY  chr6    78172677-78172685       CCDS4986.1      HTR1B   -
  chr6:g.78172677_78172685/c.436_444/p.D146_Y148  cds_in_exon_1
  protein_sequence=DRY;cDNA_sequence=GACCGCTAC;gDNA_sequence=GTAGCGGTC
```
One can also use wildcard `x` (lowercase) in the motif.
```
#!bash
transvar revanno -i 'HTR1B.p.365_369refNPxxY' --ccds
```
```
#!text
HTR1B.p.365_369refNPxxY chr6    78172014-78172028       CCDS4986.1      HTR1B   -
  chr6:g.78172014_78172028/c.1093_1107/p.N365_Y369        cds_in_exon_1
  protein_sequence=NPIIY;cDNA_sequence=AAC..TAT;gDNA_sequence=ATA..GTT
```
---
#### reverse annotation of protein range

```
#!bash
transvar revanno --ccds -i 'ABCB11:p.200_400'
```
outputs
```
#!text
ABCB11:p.200_400        chr2    169833195-169851872     CCDS46444.1     ABCB11  -
  chr2:g.169833195_169851872/c.598_1200/p.T200_K400       cds_in_exons_[6,7,8,9,10,11]
  protein_sequence=TRF..DRK;cDNA_sequence=ACA..AAA;gDNA_sequence=TTT..TGT
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
PIK3CA:E545K    chr3    178936091-178936092-178936093   ENST00000263967 PIK3CA  +
  chr3:g.178936091G>A/c.1633G>A/p.E545K   cds_in_exon_10
  reference_codon=GAG;candidate_codons=AAG,AAA;
  candidate_mnv_variants=chr3:g.178936091_178936093GAG>AAA;missense
```

One may encounter **ambiguous cases** where the multiple substitutions exist in explaining the amino acid change. For example,
```
#!bash
transvar revanno -i ACSL4:p.R133R --ccds
```
```
#!text
ACSL4:p.R133R   chrX    108926078-108926080     CCDS14548.1 (protein_coding)    ACSL4   -
  chrX:g.108926078G>T/c.399C>A/p.R133R
  cds_in_exon_2
  reference_codon=CGC;candidate_codons=AGG,AGA,CGA,CGG,CGT;
  candidate_snv_variants=chrX:g.108926078G>C,chrX:g.108926078G>A;
  candidate_mnv_variants=chrX:g.108926078_108926080GCG>CCT,chrX:g.108926078_108926080GCG>TCT;
  synonymous
```
In those cases, TransVar prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

If one wishes, he/she can also turn on the `--dbsnp` option to show further potential match to the dbSNP database. (This requires 1) pysam library installed in python path; 2) dbsnp library downloaded via `transvar config --download_hg19_dbsnp`)
```
#!bash
transvar revanno -i 'A1CF:p.A309A' --ccds --dbsnp
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
PIK3CA:c.1633G>A        chr3    178936091-178936093     CCDS43171.1 (protein_coding)    PIK3CA  +
  chr3:g.178936091G>A/c.1633G>A/p.E545K   cds_in_exon_9
  missense;reference_codon=GAG;alternative_codon=AAG
```

The SNV can be in the intronic region, e.g.,
```
#!bash
transvar revanno --ccds -i 'ABCB11:c.1198-8C>A'
```
outputs
```
#!text
ABCB11:c.1198-8C>A      chr2    169833205       CCDS46444.1 (protein_coding)    ABCB11  -
  chr2:g.169833205G>T/c.1198-8C>A/.       intron_between_exon_10_and_11   .
```
---

#### reverse annotation of cDNA region

```
#!bash
transvar revanno --ccds -i 'ABCB11:c.1198-8_1202'
```
outputs
```
#!text
ABCB11:c.1198-8_1202    chr2    169833193-169833205     CCDS46444.1     ABCB11  -
  chr2:g.169833193_169833205GGTTTCTGGAGTG/c.1198-8_1202CACTCCAGAAACC/p.400_401KP
  from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
  acceptor_splice_site_on_exon_12_included
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
ACIN1:c.1932_1933insATTCAC      chr14   23548785        CCDS9587.1 (protein_coding)     ACIN1   -
  chr14:g.23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH
  inside_[cds_in_exon_6]
  left_align_gDNA=g.23548785_23548786insGTGAAT;unalign_gDNA=g.23548785_23548786insGTGAAT;
  insertion_gDNA=GTGAAT;
  left_align_cDNA=c.1932_1933insATTCAC;unalign_cDNA=c.1932_1933insATTCAC;
  insertion_cDNA=ATTCAC;
  left_align_protein=p.R644_S645insIH;unalign_protein=p.R644_S645insIH;
  phase=0
```
`Phase = 0,1,2` indicates whether the insertion happen after the 3rd, 1st or 2nd base of a codon, respectively. An insertion *in phase* refers to one with `Phase=0`.

Example: to annotate an **out-of-phase, in-frame insertion**,
```
#!bash
transvar revanno --ccds -i 'ACIN1:c.1930_1931insATTCAC'
```
```
#!text
ACIN1:c.1930_1931insATTCAC      chr14   23548787        CCDS9587.1 (protein_coding)     ACIN1   -
  chr14:g.23548792_23548793insTGTGAA/c.1930_1931insATTCAC/p.S643_R644insHS
  inside_[cds_in_exon_6]
  left_align_gDNA=g.23548787_23548788insGTGAAT;unalign_gDNA=g.23548787_23548788insGTGAAT;
  insertion_gDNA=TGTGAA;
  left_align_cDNA=c.1925_1926insTTCACA;unalign_cDNA=c.1930_1931insATTCAC;
  insertion_cDNA=ATTCAC;
  left_align_protein=p.R642_S643insSH;unalign_protein=p.S643_R644insHS;
  phase=1
```
Reverse annotation can result in different identifiers after left/right alignments, e.g., 
```
#!bash
transvar revanno --ccds -i 'AATK:c.3976_3977insCGCCCA'
```
results in
```
AATK:c.3976_3977insCGCCCA       chr17   79093287        CCDS45807.1 (protein_coding)    AATK    -
  chr17:g.79093282_79093287dupTGGGCG/c.3988_3993dupACGCCC/p.T1330_P1331dupTP
  inside_[cds_in_exon_13]
  left_align_gDNA=g.79093270_79093271insGGGCGT;unalign_gDNA=g.79093282_79093287dupTGGGCG;
  insertion_gDNA=TGGGCG;
  left_align_cDNA=c.3976_3977insCGCCCA;unalign_cDNA=c.3976_3977insCGCCCA;
  insertion_cDNA=ACGCCC;
  left_align_protein=p.A1326_P1327insPT;unalign_protein=p.A1326_P1327insPT;phase=1
```
Note how insertion switch to duplication when 5'flanking is identical. This conforms to HGVS recommendation to replace insertion notation with duplication when possible.

Example: to annotate a **frame-shift insertion**, frameshift mutations have not alternative alignments. Hence only cDNA and gDNA have left alignment and unalignment reports.
```
#!bash
transvar revanno --ccds -i 'AAAS:c.1225_1226insG'
```
results in
```
#!text
AAAS:c.1225_1226insG    chr12   53702089        CCDS8856.1 (protein_coding)     AAAS    -
  chr12:g.53702093dupC/c.1225dupG/p.E409Gfs*17
  inside_[cds_in_exon_13]
  left_align_gDNA=g.53702089_53702090insC;unalign_gDNA=g.53702089_53702090insC;
  insertion_gDNA=C;
  left_align_cDNA=c.1221_1222insG;unalign_cDNA=c.1225dupG;
  insertion_cDNA=G
```

Example: to annotate an **intronic insertion**,
```
#!bash
transvar revanno --ccds -i 'ADAM33:c.991-3_991-2insC'
```
outputs
```
#!text
ADAM33:c.991-3_991-2insC        chr20   3654145 CCDS13058.1 (protein_coding)    ADAM33  -
  chr20:g.3654151dupG/c.991-3dupC/.
  inside_[intron_between_exon_10_and_11]
  left_align_gDNA=g.3654145_3654146insG;unalign_gDNA=g.3654145_3654146insG;
  insertion_gDNA=G;
  left_align_cDNA=c.991-9_991-8insC;unalign_cDNA=c.991-3dupC;
  insertion_cDNA=C
```
In the case of intronic insertions, amino acid identifier is not applicable, represented in a `.`. But cDNA and gDNA identifier are right-aligned according to their natural order, respecting HGVS nomenclature.

Insertion could occur to *splice sites*. TransVar identifies such cases and report splice site and repress translation of protein change.
```
#!bash
transvar revanno --ccds -i 'ADAM33:c.991_992insC'
```
results in
```
#!text
ADAM33:c.991_992insC    chr20   3654142 CCDS13058.1 (protein_coding)    ADAM33  -
  chr20:g.3654142_3654143insG/c.991_992insC/.
  inside_[cds_in_exon_11]
  left_align_gDNA=g.3654142_3654143insG;unalign_gDNA=g.3654142_3654143insG;insertion_gDNA=G;
  left_align_cDNA=c.991_992insC;unalign_cDNA=c.991_992insC;insertion_cDNA=C;
  acceptor_splice_site_on_exon_12
```

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

A4GNT:c.694_696delTTG   chr3    137843435-137843437     CCDS3097.1 (protein_coding)     A4GNT   -
  chr3:g.137843435_137843437delACA/c.694_696delTTG/p.L232del
  inside_[cds_in_exon_2]
  left_align_gDNA=g.137843433_137843435delCAA;unaligned_gDNA=g.137843433_137843435delCAA;
  left_align_cDNA=c.692_694delTGT;unalign_cDNA=c.694_696delTTG;
  left_align_protein=p.L232del;unalign_protein=p.L232del;
  deletion_gDNA=ACA;deletion_cDNA=TTG
```

Example: to annotate a **in-frame, out-of-phase deletion**,
```
#!bash
transvar revanno --ccds -i 'ABHD15.c.431_433delGTG'
```
```
#!text
ABHD15.c.431_433delGTG  chr17   27893552-27893554       CCDS32602.1 (protein_coding)    ABHD15  -
  chr17:g.27893552_27893554delCAC/c.431_433delGTG/p.C144_V145delinsF
  inside_[cds_in_exon_1]
  left_align_gDNA=g.27893552_27893554delCAC;unaligned_gDNA=g.27893552_27893554delCAC;
  left_align_cDNA=c.431_433delGTG;unalign_cDNA=c.431_433delGTG;
  deletion_gDNA=CAC;deletion_cDNA=GTG
```

Example: to annotate a **frame-shift deletion**,
```
#!bash
transvar revanno --ccds -i 'AADACL3.c.374delG'
```
```
#!text
AADACL3.c.374delG       chr1    12785494-12785494       CCDS41252.1 (protein_coding)    AADACL3 +
  chr1:g.12785494delG/c.374delG/p.C125Ffs*17      cds_in_exon_3
  left_align_gDNA=g.12785494delG;unaligned_gDNA=g.12785494delG;
  left_align_cDNA=c.374delG;unalign_cDNA=c.374delG;
  deletion_gDNA=G;deletion_cDNA=G
```

Example: to annotate a **deletion that span from intronic to coding region**, protein prediction is suppressed due to loss of splice site.
```
#!bash
transvar revanno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'
```
```
#!text
ABCB11:c.1198-8_1199delcactccagAA       chr2    169833196-169833205
  CCDS46444.1 (protein_coding)    ABCB11  -
  chr2:g.169833196_169833205delTTCTGGAGTG/c.1198-8_1199delCACTCCAGAA/.
  from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
  left_align_gDNA=g.169833196_169833205delTTCTGGAGTG;
  unaligned_gDNA=g.169833196_169833205delTTCTGGAGTG;
  left_align_cDNA=c.1198-8_1199delCACTCCAGAA;
  unalign_cDNA=c.1198-8_1199delCACTCCAGAA;
  acceptor_splice_site_on_exon_12_lost;
  deletion_gDNA=TTCTGGAGTG;deletion_cDNA=CACTCCAGAA
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
A1CF:c.508_509CC>TT     chr10   52595929-52595930       CCDS7241.1 (protein_coding)     A1CF    -
  chr10:g.52595929_52595930GG>AA/c.508_509CC>TT/p.P170L
  inside_[cds_in_exon_4]  codon_cDNA=508-509-510
```

Block substitution does not necessarily results in block substitution in amino acid. For example, the following substitution results in a deletion, where protein alternative alignment should be reported.
```
#!bash
transvar revanno --ccds -i 'CSRNP1.c.1212_1224>GGAGGAGGAA'
```
```
#!text
CSRNP1.c.1212_1224>GGAGGAGGAA   chr3    39185092-39185104       CCDS2682.1 (protein_coding)
  CSRNP1  -
  chr3:g.39185092_39185104TTCCTCCTCCTCC>TTCCTCCTCC/c.1212_1224GGAGGAGGAGGAA>GGAGGAGGAA/p.E411del
  inside_[cds_in_exon_4]
  begin_codon_cDNA=1210-1211-1212;end_codon_cDNA=1222-1223-1224;
  left_align_protein=p.E405del;unalign_protein=p.E408del
```

Likewise, block substitution could occur to **intronic region**,
```
#!bash
transvar revanno --ccds -i 'A1CF:c.1460+2_1460+3TG>CC'
```
```
#!text
A1CF:c.1460+2_1460+3TG>CC       chr10   52570797-52570798       CCDS7241.1 (protein_coding)
  A1CF    -       chr10:g.52570797_52570798CA>GG/c.1460+2_1460+3TG>CC/.
  inside_[intron_between_exon_9_and_10]   .
```

When block substitution occurs **across splice site**, TransVar put a tag in the info fields and does not predict amino acid change.
```
#!bash
transvar revanno --ccds -i 'A1CF:c.1459_1460+3ATGTG>CC'
```
```
#!text

A1CF:c.1459_1460+3ATGTG>CC      chr10   52570797-52570801       CCDS7241.1 (protein_coding)  A1CF    -
  chr10:g.52570797_52570801CACAT>GG/c.1459_1460+3ATGTG>CC/.
  from_[intron_between_exon_9_and_10]_to_[cds_in_exon_9]
  donor_splice_site_on_exon_10
```

---

#### reverse annotation of nucleotide duplication

Duplication can be thought of as special insertion where the inserted sequence is identical to the sequence flanking the breakpoint.
Similar to insertion, the annotation of duplication may possess alternative alignment.

Example: to annotate a duplication coding region,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1669_1674dup'
```
```
#!text

CHD7:c.1669_1674dup     chr8    61693569        CCDS47865.1 (protein_coding)    CHD7    +
  chr8:g.61693564_61693569dupCCCGTC/c.1669_1674dup/p.P558_S559dupPS
  inside_[cds_in_exon_2]
  left_align_gDNA=g.61693561_61693562insTCCCCG;unalign_gDNA=g.61693562_61693567dupTCCCCG;
  insertion_gDNA=CCCGTC;
  left_align_cDNA=c.1668_1669insTCCCCG;unalign_cDNA=c.1669_1674dupTCCCCG;
  insertion_cDNA=CCCGTC;
  left_align_protein=p.H556_S557insSP;unalign_protein=p.S557_P558dupSP;
  phase=0
```

Example: a duplication on the nucleotide level may lead to frame-shift or block substitution on the amino acid level,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1668_1669dup'
```
```
#!text
CHD7:c.1668_1669dup     chr8    61693562        CCDS47865.1 (protein_coding)    CHD7    +
  chr8:g.61693561_61693562dupTT/c.1668_1669dup/p.S557Ffs*8
  inside_[cds_in_exon_2]
  left_align_gDNA=g.61693560_61693561insTT;unalign_gDNA=g.61693561_61693562dupTT;
  insertion_gDNA=TT;
  left_align_cDNA=c.1667_1668insTT;unalign_cDNA=c.1668_1669dupTT;
  insertion_cDNA=TT
```

Example: to annotate a duplication in intronic region,
```
#!bash
transvar revanno --ccds -i 'CHD7:c.1666-5_1666-3dup'
```
```
#!text
CHD7:c.1666-5_1666-3dup chr8    61693556        CCDS47865.1 (protein_coding)    CHD7    +
  chr8:g.61693554_61693556dupCTC/c.1666-5_1666-3dup/.
  inside_[intron_between_exon_1_and_2]
  left_align_gDNA=g.61693553_61693554insCTC;unalign_gDNA=g.61693554_61693556dupCTC;
  insertion_gDNA=CTC;
  left_align_cDNA=c.1666-6_1666-5insCTC;unalign_cDNA=c.1666-5_1666-3dupCTC;
  insertion_cDNA=CTC
```

---

#### reverse annotation of amino acid insertion

```
#!bash
transvar revanno --ccds -i 'AATK.p.P1331_A1332insTP'
```
```
#!text

AATK.p.P1331_A1332insTP chr17   79093267        CCDS45807.1 (protein_coding)    AATK    -
  chr17:g.(79093267ins6)/c.(3997_3991ins6)/p.T1330_P1331dupTP
  cds_in_exon_13
  left_align_protein=p.A1326_P1327insPT;unalign_protein=p.T1330_P1331dupTP;
  insertion_cDNA=ACACCT;insertion_gDNA=AGGTGT;imprecise
```

#### reverse annotation of amino acid deletion
```
#!bash
transvar revanno --ccds -i 'AADACL4.p.W263_I267delWRDAI'
```
```
#!text
AADACL4.p.W263_I267delWRDAI     chr1    12726309-12726323       CCDS30590.1 (protein_coding)    AADACL4 +
  chr1:g.12726309_12726323/c.787_801/p.W263_I267delWRDAI
  inside_[cds_in_exon_4]
  left_align_protein=p.W263_I267delWRDAI;unalign_protein=p.W263_I267delWRDAI;imprecise
```

#### reverse annotation of amino acid block substitution
```
#!bash
transvar revanno --ccds -i 'ABCC3:p.Y556_V557delinsRRR'
```
```
#!text
ABCC3:p.Y556_V557delinsRRR      chr17   48745254-48745259       CCDS32681.1 (protein_coding)   ABCC3   +
  chr17:g.48745254_48745259TACGTG>AGGAGGAGG/c.1666_1671TACGTG>AGGAGGAGG/p.Y556_V557delinsRRR
  cds_in_exon_13  imprecise
```

#### reverse annotation of amino acid frame-shift

```
#!bash
transvar revanno --ccds -i 'A1BG.p.G132fs*2'
```
```
#!text
A1BG.p.G132fs*2 chr19   58863866-58863867-58863868      CCDS12976.1 (protein_coding)    A1BG    -
  chr19:g.58863860-58863868/c.394-402/p.G132fs*2  cds_in_exon_4   imprecise
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
DHODH:G152R  DHODH.p.G16   chr16 72050942-72050943-72050944   72050942-72050943-72050944   NM_001361.4[RefSeq]/XM_005255828.1[RefSeq]
DHODH:G152R  DHODH.p.G9    chr16 72050942-72050943-72050944   72050942-72050943-72050944   NM_001361.4[RefSeq]/XM_005255829.1[RefSeq]
DHODH:G152R  DHODH.p.G124  chr16 72050942-72050943-72050944   72050942-72050943-72050944   NM_001361.4[RefSeq]/XM_005255827.1[RefSeq]
```
TransVar outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

Example: to run `transvar codonsearch` to **batch process** a list of mutation identifiers.
```
#!bash
transvar codonsearch -l example/input_table2 --ccds -m 1 -o 1
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
CDKN2A.p61  CDKN2A.p.76   chr9  21971175-21971176-21971177  21971173-21971174-21971175    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p61  CDKN2A.p.75   chr9  21971175-21971176-21971177  21971176-21971177-21971178    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69  CDKN2A.p.54   chr9  21971194-21971195-21971196  21971196-21971197-21971198    CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69  CDKN2A.p.55   chr9  21971194-21971195-21971196  21971193-21971194-21971195    CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69  CDKN2A.p.83   chr9  21971151-21971152-21971153  21971152-21971153-21971154    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69  CDKN2A.p.84   chr9  21971151-21971152-21971153  21971149-21971150-21971151    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69  CDKN2A.p.54   chr9  21971194-21971195-21971196  21971196-21971197-21971198    CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69  CDKN2A.p.55   chr9  21971194-21971195-21971196  21971193-21971194-21971195    CCDS6511.2[CDDS]/CCDS6510.1[CDDS],CCDS6511.2[CDDS]/CCDS56565.1[CDDS]
CDKN2A.p69  CDKN2A.p.83   chr9  21971151-21971152-21971153  21971152-21971153-21971154    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
CDKN2A.p69  CDKN2A.p.84   chr9  21971151-21971152-21971153  21971149-21971150-21971151    CCDS6510.1[CDDS]/CCDS6511.2[CDDS],CCDS56565.1[CDDS]/CCDS6511.2[CDDS]
ERBB2.p755  ERBB2.p.785   chr17 37881024-37881025-37881026  37881024-37881025-37881026    CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
ERBB2.p755  ERBB2.p.725   chr17 37880219-37880220-37880221  37880219-37880220-37880221    CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
ERBB2.p755  ERBB2.p.785   chr17 37881024-37881025-37881026  37881024-37881025-37881026    CCDS45667.1[CDDS]/CCDS32642.1[CDDS]
ERBB2.p755  ERBB2.p.725   chr17 37880219-37880220-37880221  37880219-37880220-37880221    CCDS32642.1[CDDS]/CCDS45667.1[CDDS]
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
MET.p.1010  MET.p.1047    7   116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/XM_005250353.1[RefSeq]
MET.p.1010  MET.p.580     7   116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/XM_005250354.1[RefSeq]
MET.p.1010  MET.p.991     7   116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/NM_001127500.1[RefSeq]
MET.p.1010  MET.p.1029    7   116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/XM_005250353.1[RefSeq]
MET.p.1010  MET.p.543     7   116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/XM_005250354.1[RefSeq]
MET.p.1010  MET.p.992     7   116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/NM_000245.2[RefSeq]
MET.p.1010  MET.p.973     7   116411932-116411933-116411934   116411932-116411933-116411934   XM_005250353.1[RefSeq]/NM_000245.2[RefSeq]
MET.p.1010  MET.p.1028    7   116412043-116414935-116414936   116412043-116414935-116414936   NM_000245.2[RefSeq]/NM_001127500.1[RefSeq]
MET.p.1010  MET.p.562     7   116411989-116411990-116411991   116411989-116411990-116411991   NM_001127500.1[RefSeq]/XM_005250354.1[RefSeq]
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
chr3:178936091.G>A      chr3    178936091       CCDS43171.1 (protein_coding)    PIK3CA  +
  chr3:g.178936091G>A/c.1633G>A/p.E545K
  cds_in_exon_9
  missense;codon_pos=178936091-178936092-178936093;ref_codon_seq=GAG
```

Another example:
```
#!bash
transvar anno -i "9:135782704C>G" --ccds
```
outputs
```
#!text
9:135782704C>G  chr9    135782704       CCDS6956.1 (protein_coding)     TSC1    -
  chr9:g.135782704C>G/c.1317G>C/p.L439L   cds_in_exon_11
  synonymous;codon_pos=135782704-135782705-135782706;ref_codon_seq=CTG
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
chr3:g.178936091_178936192      chr3    178936091-178936192     CCDS43171.1 (protein_coding)    PIK3CA  +
  chr3:g.178936091_178936192/c.1633_1664+70/p.E545_R555
  from_[cds_in_exon_9]_to_[intron_between_exon_9_and_10]
  donor_splice_site_on_exon_8_included;
  start_codon=178936091-178936092-178936093;end_codon=178936121-178936122-178936984
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
9:g.133750356_137990357 chr9    133750356-137990357     CCDS35165.1 (protein_coding),CCDS6986.1 (protein_coding) . .
  chr9:g.133750356_137990357/./.
  from_[cds_in_exon_7_(ABL1)]_to_[intron_between_exon_4_and_5_(OLFM1)]_spanning_[51_genes] .
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
9:g.133750356_1337503570    chr9    133750356-141213431     CCDS35165.1 (protein_coding),   .       .
  chr9:g.133750356_141213431/./.
  from_[cds_in_exon_7;ABL1]_to_[intergenic_between_EHMT1(484,026_bp_downstream)_and_3'-telomere(0_bp)]_spanning_[136_genes]       .
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
chr2:234183368_234183380del     chr2    234183368-234183380     CCDS2502.2 (protein_coding)     ATG16L1 +
  chr2:g.234183368_234183380del13/c.841_853del13/p.T281Lfs*5      inside_[cds_in_exon_8]
  left_align_gDNA=g.234183367_234183379del13;unaligned_gDNA=g.234183368_234183380del13;
  left_align_cDNA=c.840_852del13;unalign_cDNA=c.841_853del13
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
chr2:234183368_234183379del     chr2    234183368-234183379     CCDS2502.2 (protein_coding)     ATG16L1 +
  chr2:g.234183368_234183379del12/c.841_852del12/p.T281_G284delTHPG
  inside_[cds_in_exon_8]
  left_align_gDNA=g.234183367_234183378del12;unaligned_gDNA=g.234183368_234183379del12;
  left_align_cDNA=c.840_851del12;unalign_cDNA=c.841_852del12;
  left_align_protein=p.T281_G284delTHPG;unalign_protein=p.T281_G284delTHPG
```

Another example
```
#!bash
transvar anno --ccds -i '12:53703425_53703427del'
```
outputs
```
#!text
12:53703425_53703427del chr12   53703427-53703429       CCDS8856.1 (protein_coding)     AAAS    -
  chr12:g.53703427_53703429delCCC/c.769_771delGGG/p.G257del
  inside_[cds_in_exon_8]
  left_align_gDNA=g.53703424_53703426delCCC;unaligned_gDNA=g.53703425_53703427delCCC;
  left_align_cDNA=c.766_768delGGG;unalign_cDNA=c.768_770delGGG;
  left_align_protein=p.G256del;unalign_protein=p.G256del
```
Note the difference between left and right-aligned identifiers on both protein level and cDNA level.

An in-frame out-of-phase deletion
```
#!bash
transvar anno -i "chr2:234183372_234183383del" --ccds
```
outputs
```
chr2:234183372_234183383del     chr2    234183372-234183383     CCDS2502.2 (protein_coding)     ATG16L1 +
  chr2:g.234183372_234183383del12/c.845_856del12/p.H282_G286delinsR
  inside_[cds_in_exon_8]
  left_align_gDNA=g.234183372_234183383del12;unaligned_gDNA=g.234183372_234183383del12;
  left_align_cDNA=c.845_856del12;unalign_cDNA=c.845_856del12
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

2:69741762insTGC        chr2    69741762        CCDS1893.2 (protein_coding)     AAK1    -
  chr2:g.69741780_69741782dupCTG/c.1614_1616dupGCA/p.Q546dupQ
  cds_in_exon_12
  left_align_gDNA=g.69741762_69741763insTGC;unalign_gDNA=g.69741762_69741763insTGC;
  insertion_gDNA=CTG;
  left_align_cDNA=c.1596_1597insCAG;unalign_cDNA=c.1614_1616dupGCA;
  insertion_cDNA=GCA;
  left_align_protein=p.Y532_Q533insQ;unalign_protein=p.Q539dupQ;
  phase=2
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
7:121753754insCA        chr7    121753754       CCDS5783.1 (protein_coding)     AASS    -
  chr7:g.121753754_121753755insCA/c.1064_1065insGT/p.I355Mfs*10
  cds_in_exon_9
  left_align_gDNA=g.121753753_121753754insAC;unalign_gDNA=g.121753754_121753755insCA;
  insertion_gDNA=CA;
  left_align_cDNA=c.1063_1064insTG;unalign_cDNA=c.1063_1064insTG;
  insertion_cDNA=GT
```

```
#!bash
transvar anno -i '17:79093270insGGGCGT' --ccds
```
outputs
```
#!text
17:79093270insGGGCGT    chr17   79093270        CCDS45807.1 (protein_coding)    AATK    -
  chr17:g.79093282_79093287dupTGGGCG/c.3988_3993dupACGCCC/p.T1330_P1331dupTP
  cds_in_exon_13
  left_align_gDNA=g.79093270_79093271insGGGCGT;unalign_gDNA=g.79093270_79093271insGGGCGT;
  insertion_gDNA=TGGGCG;
  left_align_cDNA=c.3976_3977insCGCCCA;unalign_cDNA=c.3988_3993dupACGCCC;
  insertion_cDNA=ACGCCC;
  left_align_protein=p.A1326_P1327insPT;unalign_protein=p.T1330_P1331dupTP;
  phase=0
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
7:117230474insA chr7    117230474       CCDS5773.1 (protein_coding)     CFTR    +
  chr7:g.117230479dupA/c.1752dupA/p.E585Rfs*4
  cds_in_exon_13
  left_align_gDNA=g.117230474_117230475insA;unalign_gDNA=g.117230474_117230475insA;
  insertion_gDNA=A;
  left_align_cDNA=c.1747_1748insA;unalign_cDNA=c.1747_1748insA;
  insertion_cDNA=A
```
Notice the right alignment of cDNA level insertion and the left alignment reported as additional information.

A in-frame, in-phase insertion
```
#!bash
transvar anno -i '12:109702119insACC' --ccds
```
```
#!text
12:109702119insACC      chr12   109702119       CCDS31898.1 (protein_coding)    ACACB   +
  chr12:g.109702119_109702120insACC/c.6870_6871insACC/p.Y2290_H2291insT
  cds_in_exon_49
  left_align_gDNA=g.109702118_109702119insCAC;unalign_gDNA=g.109702119_109702120insACC;
  insertion_gDNA=ACC;
  left_align_cDNA=c.6869_6870insCAC;unalign_cDNA=c.6870_6871insACC;
  insertion_cDNA=ACC;
  left_align_protein=p.Y2290_H2291insT;unalign_protein=p.Y2290_H2291insT;
  phase=0
```

#### annotate block substitution from genomic locations

A block-substitution that results in a frameshift.
```
#!bash
transvar anno -i 'chr10:g.27329002_27329002A>AT' --ccds
```
```
#!text
chr10:g.27329002_27329002A>AT   chr10   27329002-27329002       CCDS41499.1 (protein_coding)    ANKRD26 -
  chr10:g.27329002_27329002A>AT/c.2267_2267T>AT/p.M756Nfs*6       cds_in_exon_21  .
```

A block-substitution that is in-frame,
```
#!bash
transvar anno -i 'chr10:g.52595929_52595930GG>AA' --ccds
```
```
#!text
chr10:g.52595929_52595930GG>AA  chr10   52595929-52595930       CCDS7241.1 (protein_coding)     A1CF    -
  chr10:g.52595929_52595930GG>AA/c.508_509CC>TT/p.P170L   inside_[cds_in_exon_4]  codon_cDNA=508-509-510
```

### annotate promoter region

One can define the promoter boundary through the `--prombeg` and `--promend` option. Default promoter region is defined from 1000bp upstream of the transcription start site to the transcription start site. One could customize this setting to e.g., [-1000bp, 2000bp] by

```
#!bash
transvar anno -i 'chr19:41950335_41951908' --ensembl --prombeg 2000 --promend 1000 --refversion mm10
```
```
#!text
chr19:41950335_41951908 chr19   41950335-41951908       ENSMUST00000166517      MMS19   -
chr19:g.41950335_41951908/c.1-564_593+417/p.L1_E198
from_[intron_between_exon_1_and_2]_to_[intergenic_between_MMS19(564_bp_upstream)_and_MMS19(1,189_bp_downstream)]
promoter_overlap_1564_bp(99.43%);start_codon=41950753-41950752-41950083;end_codon=41951344-41951343-41951342
```
The result shows that 99.43% of the target region is inside the promoter region. The overlap is as long as 1564 base pairs.

#### annotate non-coding RNA
Given Ensembl, GENCODE or RefSeq database, one could annotate non-coding transcripts such as lncRNA.
E.g.,
```
#!bash
transvar anno --gencode -i 'chr1:3985200_3985300' -refversion mm10
```
results in
```
#!text
chr1:3985200_3985300    chr1    3985200-3985300 ENSMUST00000194643.1 (lincRNA)
  RP23-333I7.1    -       chr1:g.3985200_3985300/c.121_221/.      inside_[noncoding_exon_2]       .
```
or
```
#!bash
transvar anno --refseq -i 'chr14:20568338_20569581' -refversion mm10
```
results in
```
#!text
chr14:20568338_20569581 chr14   20568338-20569581       NR_033571.1 (lncRNA)
  1810062O18Rik   +       chr14:g.20568338_20569581/c.260-1532_260-289/.  inside_[intron_between_exon_4_and_5]       .
```

or using Ensembl
```
#!bash
transvar anno --ensembl -i 'chr1:29560_29570'
```
results in
```
#!text
chr1:29560_29570        chr1    29560-29570     ENST00000473358 (lincRNA)       MIR1302-10      +
  chr1:g.29560_29570/c.7_17/.     inside_[noncoding_exon_1]       .
```

### FAQ

+ I got 'GeneNotRecognized', what's wrong?

Most likely you forgot to specify a transcipt definition such as `--ccds` or `--ensembl`. Sometimes there are non-canonical names for genes, this can be fixed through the `--alias` option and specify an alias table. TransVar comes with alias table from UCSC knownGene.

### Technical notes

TransVar follows in full the HGVS nomenclature while annotating protein level mutation identifiers. For example, a out-of-phase, in frame insertion, `ACIN1:c.1930_1931insATTCAC` will be annotated with `p.S643_R644insHS` rather than `R644delinsHSR`. Protein level mutation will be generated as if no nucleotide mutation information exists.

## Future work

 + add cytoband annotation
 + option to output full long deletion sequence
 + imprecise annotation
 + forward annotation of binding sites
 + forward annotation of structural variation breakpoints
 + begin codon and end codon in deletion
 + distinguish non-transcribable element and suppress promoter setting (like "retained intron")
 + dbsnp id 

## Bug report and feature request

If you find any bug (very likely due to the complexity of genomics:-)) or you wish any feature, please direct to Wanding Zhou <zhouwanding@gmail.com>. Thank you.

## Reference

We are working on an application note on this topic :-).

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.


