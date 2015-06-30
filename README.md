**TransVar** is a versatile annotator for 3-way conversion and annotation among genomic characterization(s) of mutations (e.g., `chr3:g.178936091G>A`) and transcript-dependent annotation(s) (e.g., `PIK3CA:p.E545K` or `PIK3CA:c.1633G>A`, or `NM_006218.2:p.E545K`). It is particularly designed with the functionality of resolving ambiguous mutation annotations arising from differential transcript usage. TransVar keeps awareness of the underlying unknown transcript structure (exon boundary, reference amino acid/base) while performing reverse annotation (from protein level to cDNA level).
TransVar has the following features:

 + supports HGVS nomenclature
 + supports input from gene name, transcript ID, UniProt ID and other aliases.
 + supports both left-alignment and right-alignment convention in reporting indels and duplications.
 + supports annotation of a region based on a transcript-dependent characterization
 + supports mutations at both coding region and intronic/UTR regions
 + supports noncoding RNA annotation
 + supports VCF inputs
 + supports long haplotype decomposition
 + supports single nucleotide variation (SNV), insertions and deletions (indels) and block substitutions
 + supports transcript annotation from commonly-used databases such as Ensembl, NCBI RefSeq and GENCODE etc
 + supports GRCh36, 37, 38 (human),  GRCm38 (mouse), NCBIM37 (mouse)
 + supports >60 other genomes available from Ensembl
 + functionality of forward annotation.

--------

[TOC]

--------

### Download and install

#### dependency

requires just Python 2.7. 

#### download the program

current stable version: [v2.1.2.20150630](https://bitbucket.org/wanding/transvar/get/v2.1.2.20150630.zip)

For previous versions, see [TAGS](https://bitbucket.org/wanding/transvar/overview#tags).

for also stable 2.0.x version:
[v2.0.12.20150626](https://bitbucket.org/wanding/transvar/get/v2.0.12.20150626.zip)

for older 1.x version:
[v1.40](https://bitbucket.org/wanding/transvar/get/v1.40.zip)

#### install

##### System-wise install (need root)
```
#!bash
sudo python setup.py install
```

##### Local install
```
#!bash
python setup.py install --prefix [folder]
```
After install, there will be two subfolders in `[folder]/lib` (which would contain libraries) and `[folder]/bin` (which would contain transvar executable).
When you run transvar, make sure `[folder]/lib/python2.7/site-packages` is in your PYTHONPATH. In some occasions, you need to `mkdir -p [folder]/lib/python2.7/site-packages` to make sure it exists before you could run `setup.py`.
You can add it by putting
`export PYTHONPATH=$PYTHONPATH:[folder]/lib/python-2.7/site-packages/` to your `.bashrc` (or `.profile` depending on your OS).

The installed executable is `[folder]/bin/transvar`.

#### quick start

Here we show how one can use TransVar on human hg19 (GRCh37).
```
#!bash
# set up databases
transvar config --download_anno --refversion hg19
# in case, if you don't have a reference
transvar config --download_ref --refversion hg19
# in case you do have a reference to link
transvar config -k reference -v [path_to_hg19.fa] --refversion hg19

$ transvar panno -i 'PIK3CA:p.E545K' --ucsc --ccds
```
outputs show two hits from the two databases, i.e., UCSC and CCDS.
```
#!text
PIK3CA:p.E545K	NM_006218 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_10
   reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_variants=chr3:g.17
   8936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G>A);missense
PIK3CA:p.E545K	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_9
   reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_variants=chr3:g.17
   8936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G>A);missense
```
One could provide input based on transcript ID, e.g `NM_006218.1:p.E545K` and TransVar would automatically restrict to the provided transcript.
```
$ transvar panno -i 'NM_006218.2:p.E545K' --ucsc --ccds
```
outputs
```
#!text
NM_006218.2:p.E545K	NM_006218 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_10
   reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_variants=chr3:g.17
   8936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G>A);missense
```

#### install/specify reference genome assembly

For some genome assembly we provide download via `transvar config --download_ref --refversion hg19`. This will download the faidx indexed reference.
For other genome assemblies,one could download the genome and index it by, e.g.,
```
#!bash
transvar index --reference hg19.fa
```
Under the hood, TransVar uses the `samtools faidx`. So one could use any existing faidx indices without a glitch.
Once downloaded and indexed, the genome can be used through the "--reference" option followed by path to the genome.

To set the default location of genome file for a reference version, say, to ./hg19.fa,
```
#!bash
transvar config -k reference -v ./hg19.fa --refversion hg19
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
transvar config --download_anno --refversion hg19
```
will automatically download annotation from Ensembl, RefSeq etc. to `[installdir]/lib/transvar/transvar.download` directory or your local `~/.transvar.download` if the installation directory is inaccessible.
See `transvar config -h` for downloading more versions.
These will also create default mappings under the corresponding reference version section of transvar.cfg like
```
#!text
[hg19]
ucsc = /home/wzhou1/download/hg19.ucsc.txt.gz
```

One also has the option of downloading from Ensembl collection.
```
#!bash
transvar config --download_ensembl --refversion mus_musculus
```
Without specifying the refversion, user will be prompted a collection of options to choose from.

### Usage

#### specify transcript annotation

To index transcript annotation,
```
#!bash
transvar index --ensembl GRCh37.75.gtf.gz
```
creates `GRCh37.75.gtf.gz.transvardb` together with indices files `GRCh37.75.gtf.gz.transvardb.loc_idx`, `GRCh37.75.gtf.gz.transvardb.gene_idx`, `GRCh37.75.gtf.gz.transvardb.trxn_idx` etc.

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

If one download transcripts through "transvar config", TransVar would use the downloaded definition automatically (by setting the default configuration file). For example, "--ccds" would look for the downloaded CCDS definition. One can specify non-default annotation by appending a path to the option ("--ccds CCDS.current.txt.transvardb").
To set the default annotation of a particular reference version,
```
#!bash
transvar config -k ccds -v CCDS.current.txt.transvardb --refversion hg19
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
transvar config -k ccds -v ccds.fly.txt --refversion drosophila_melanogaster
```
Will create in transvar.cfg a section like
```
#!text
[drosophila_melanogaster]
ccds = ccds.fly.txt
```

To switch to a version on the fly, one could use the "--refversion" option, e.g.,
```
#!bash
transvar panno -i 'PIK3CA:p.E545K' --ucsc --refversion hg38
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
specifying `--refversion` displays the information under that reference version (without changing the default reference version setup).

#### batch processing

For all mutation types, one can batch process a list of mutation identifiers with optional transcript id to constraint the search. Take SNV for example,
```
#!bash
$ transvar panno -l example/input_table -g 1 -m 5 -t 2 --ensembl -o 2,3,4
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
ENST00000286744|15:84442328|c.243G>A	ENST00000286744 (protein_coding)	ADAMTSL3	+
   chr15:g.84442327G>A/c.242G>A/p.W81*	cds_in_exon_4
   reference_codon=TGG;candidate_codons=TAA,TAG,TGA;candidate_snv_variants=chr15
   :g.84442328G>A;candidate_mnv_variants=chr15:g.84442327_84442328delGGinsAA;mis
   sense
ENST00000286744|15:84442326|c.241T>C	ENST00000286744 (protein_coding)	ADAMTSL3	+
   chr15:g.84442326T>A/c.241T>A/p.W81R	cds_in_exon_4
   reference_codon=TGG;candidate_codons=AGG,AGA,CGA,CGC,CGG,CGT;candidate_snv_va
   riants=chr15:g.84442326T>C;candidate_mnv_variants=chr15:g.84442326_84442328de
   lTGGinsAGA,chr15:g.84442326_84442328delTGGinsCGA,chr15:g.84442326_84442328del
   TGGinsCGC,chr15:g.84442326_84442328delTGGinsCGT;missense
ENST00000369038|1:150530513|c.2270G>A	ENST00000369038 (protein_coding)	ADAMTSL4	+
   chr1:g.150530513G>A/c.2270G>A/p.G757D	cds_in_exon_12
   reference_codon=GGT;candidate_codons=GAC,GAT;candidate_mnv_variants=chr1:g.15
   0530513_150530514delGTinsAC;missense
ENST00000338316|5:7802364|c.2662G>A	ENST00000338316 (protein_coding)	ADCY2	+
   chr5:g.7802364G>A/c.2662G>A/p.V888I	cds_in_exon_21
   reference_codon=GTC;candidate_codons=ATC,ATA,ATT;candidate_mnv_variants=chr5:
   g.7802364_7802366delGTCinsATA,chr5:g.7802364_7802366delGTCinsATT;missense
ENST00000338316|5:7802365|c.2663T>C	ENST00000338316 (protein_coding)	ADCY2	+
   chr5:g.7802365T>C/c.2663T>C/p.V888A	cds_in_exon_21
   reference_codon=GTC;candidate_codons=GCA,GCC,GCG,GCT;candidate_mnv_variants=c
   hr5:g.7802365_7802366delTCinsCA,chr5:g.7802365_7802366delTCinsCG,chr5:g.78023
   65_7802366delTCinsCT;missense
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
$ transvar panno --ccds -i 'Q5VUM1:47' --uniprot
```
```
#!text
Q5VUM1:47	CCDS4972 (protein_coding)	C6ORF57	+
   chr6:g.71289191_71289193/c.139_141/p.47S	cds_in_exon_2
   protein_sequence=S;cDNA_sequence=TCC;gDNA_sequence=TCC
```
TransVar use a keyword extension `ref` in `Q5VUM1:p.47refS` to differentiate from the synonymous mutation `Q5VUM1:p.47S`. The former notation specifies that the reference protein sequence is `S` while the later specifies the target protein sequence is `S`.

---

#### reverse annotation of protein motif

For example, one can find the genomic location of a DRY motif in protein P28222 by issuing the following command,
```
#!bash
$ transvar panno -i 'P28222:p.146_148refDRY' --uniprot --ccds
```
```
#!text
P28222:p.146_148refDRY	CCDS4986 (protein_coding)	HTR1B	-
   chr6:g.78172677_78172685/c.436_444/p.D146_Y148	cds_in_exon_1
   protein_sequence=DRY;cDNA_sequence=GACCGCTAC;gDNA_sequence=GTAGCGGTC
```
One can also use wildcard `x` (lowercase) in the motif.
```
#!bash
$ transvar panno -i 'HTR1B:p.365_369refNPxxY' --ccds
```
```
#!text
HTR1B:p.365_369refNPxxY	CCDS4986 (protein_coding)	HTR1B	-
   chr6:g.78172014_78172028/c.1093_1107/p.N365_Y369	cds_in_exon_1
   protein_sequence=NPIIY;cDNA_sequence=AAC..TAT;gDNA_sequence=ATA..GTT
```
---
#### reverse annotation of protein range

```
#!bash
$ transvar panno --ccds -i 'ABCB11:p.200_400'
```
outputs
```
#!text
ABCB11:p.200_400	CCDS46444 (protein_coding)	ABCB11	-
   chr2:g.169833195_169851872/c.598_1200/p.T200_K400	cds_in_exons_[6,7,8,9,10,11]
   protein_sequence=TRF..DRK;cDNA_sequence=ACA..AAA;gDNA_sequence=TTT..TGT
```

---
#### reverse annotation of single amino acid substitution
Mutation formats acceptable in TransVar are ```PIK3CA:p.E545K``` or without reference or alternative amino acid identity, e.g., ```PIK3CA:p.545K``` or ```PIK3CA:p.E545```. TransVar takes native HGVS format inputs and outputs. The reference amino acid is used to narrow the search scope of candidate transcripts. The alternative amino acid is used to infer nucleotide change which results in the amino acid.

```
#!bash
$ transvar panno -i PIK3CA:p.E545K --ensembl
```
outputs
```
#!text
PIK3CA:p.E545K	ENST00000263967 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_10
   reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_variants=chr3:g.17
   8936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G>A);missense
```

One may encounter **ambiguous cases** where the multiple substitutions exist in explaining the amino acid change. For example,
```
#!bash
$ transvar panno -i ACSL4:p.R133R --ccds
```
```
#!text
ACSL4:p.R133R	CCDS14548 (protein_coding)	ACSL4	-
   chrX:g.108926078G>T/c.399C>A/p.R133R	cds_in_exon_2
   reference_codon=CGC;candidate_codons=AGG,AGA,CGA,CGG,CGT;candidate_snv_varian
   ts=chrX:g.108926078G>C,chrX:g.108926078G>A;candidate_mnv_variants=chrX:g.1089
   26078_108926080delGCGinsCCT,chrX:g.108926078_108926080delGCGinsTCT;synonymous
```
In those cases, TransVar prioritizes all the candidate base changes by minimizing the edit distance between the reference codon sequence and the target codon sequence. One of the optimal base changes is arbitrarily chosen as the default and all the candidates are included in the appended `CddMuts` entry.

#### annotate with additional resources

For example, one could annotate SNP with dbSNP id by downloading the dbSNP files.
This can be done by
```
#!bash
transvar config --download_dbsnp
```
TransVar automatically download dbSNP file which correspoding to the current default reference version (as set in `transvar.cfg`). This also sets the entry in `transvar.cfg`.
With dbSNP file downloaded, TransVar automatically looks for dbSNP id when performing annotation.
```
#!bash
$ transvar panno -i 'A1CF:p.A309A' --ccds
```
```
#!text
A1CF:p.A309A	CCDS7243 (protein_coding)	A1CF	-
   chr10:g.52576004T>G/c.927A>C/p.A309A	cds_in_exon_7
   reference_codon=GCA;candidate_codons=GCC,GCG,GCT;candidate_snv_variants=chr10
   :g.52576004T>C,chr10:g.52576004T>A;dbsnp=rs201831949(chr10:52576004T>G);synon
   ymous
```
Note that in order to use dbSNP, one must download the dbSNP database through `transvar config --download_dbsnp`, or by configure the `dbsnp` slot in the configure file via `transvar config -k dbsnp -v [path to dbSNP VCF]`. Manually set path for dbSNP file must have the file tabix indexed.


---

#### reverse annotation of single nucleotide variation (SNV)

TransVar infers nucleotide mutation through ```PIK3CA:c.1633G>A```. Note that nucleotide identity follows the natural sequence, i.e., if transcript is interpreted on the reverse-complementary strand, the base at the site needs to be reverse-complemented too.
```
#!bash
$ transvar canno --ccds -i 'PIK3CA:c.1633G>A'
```
outputs
```
#!text
PIK3CA:c.1633G>A	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_9
   dbsnp=rs104886003(chr3:178936091G>A);missense;reference_codon=GAG;alternative
   _codon=AAG
```

The SNV can be in the intronic region, e.g.,
```
#!bash
$ transvar canno --ccds -i 'ABCB11:c.1198-8C>A'
```
outputs
```
#!text
ABCB11:c.1198-8C>A	CCDS46444 (protein_coding)	ABCB11	-
   chr2:g.169833205G>T/c.1198-8C>A/.	intron_between_exon_10_and_11
   .
```
---

#### reverse annotation of cDNA region

```
#!bash
$ transvar canno --ccds -i 'ABCB11:c.1198-8_1202'
```
outputs
```
#!text
ABCB11:c.1198-8_1202	CCDS46444 (protein_coding)	ABCB11	-
   chr2:g.169833193_169833205GGTTTCTGGAGTG/c.1198-8_1202CACTCCAGAAACC/p.400_401KP	from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
   acceptor_splice_site_on_exon_11_at_chr2:169833198_included
```

---

#### reverse annotation of nucleotide insertion
An insertion may result in: 1) a pure insertion of amino acids; 2) a block substitution of amino acids, when insertion occur after 1st or 2nd base in a codon; or 3) a frame-shift. Following HGVS nomenclature, TransVar labels the first different amino acid and the length of the peptide util stop codon, assuming no change in the splicing.

Example: to annotate an **in-frame, in-phase insertion**,
```
#!bash
$ transvar canno --ccds -i 'ACIN1:c.1932_1933insATTCAC'
```
```
#!text
ACIN1:c.1932_1933insATTCAC	CCDS9587 (protein_coding)	ACIN1	-
   chr14:g.23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548785_23548786insGTGAAT;unalign_gDNA=g.23548785_23548786
   insGTGAAT;insertion_gDNA=GTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
   cDNA=c.1932_1933insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.R644_S64
   5insIH;unalign_protein=p.R644_S645insIH;phase=0
ACIN1:c.1932_1933insATTCAC	CCDS53889 (protein_coding)	ACIN1	-
   chr14:g.23548157_23548158insGTGAAT/c.1932_1933insATTCAC/p.P644_V645insIH	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548157_23548158insGTGAAT;unalign_gDNA=g.23548157_23548158
   insGTGAAT;insertion_gDNA=GTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
   cDNA=c.1932_1933insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.P644_V64
   5insIH;unalign_protein=p.P644_V645insIH;phase=0
ACIN1:c.1932_1933insATTCAC	CCDS55905 (protein_coding)	ACIN1	-
   chr14:g.23548785_23548786insGTGAAT/c.1932_1933insATTCAC/p.R644_S645insIH	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548785_23548786insGTGAAT;unalign_gDNA=g.23548785_23548786
   insGTGAAT;insertion_gDNA=GTGAAT;left_align_cDNA=c.1932_1933insATTCAC;unalign_
   cDNA=c.1932_1933insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.R644_S64
   5insIH;unalign_protein=p.R644_S645insIH;phase=0
```
`Phase = 0,1,2` indicates whether the insertion happen after the 3rd, 1st or 2nd base of a codon, respectively. An insertion *in phase* refers to one with `Phase=0`.

Example: to annotate an **out-of-phase, in-frame insertion**,
```
#!bash
$ transvar canno --ccds -i 'ACIN1:c.1930_1931insATTCAC'
```
```
#!text
ACIN1:c.1930_1931insATTCAC	CCDS9587 (protein_coding)	ACIN1	-
   chr14:g.23548792_23548793insTGTGAA/c.1930_1931insATTCAC/p.S643_R644insHS	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548787_23548788insGTGAAT;unalign_gDNA=g.23548787_23548788
   insGTGAAT;insertion_gDNA=TGTGAA;left_align_cDNA=c.1925_1926insTTCACA;unalign_
   cDNA=c.1930_1931insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.R642_S64
   3insSH;unalign_protein=p.S643_R644insHS;phase=1
ACIN1:c.1930_1931insATTCAC	CCDS53889 (protein_coding)	ACIN1	-
   chr14:g.23548162_23548163insAATGTG/c.1930_1931insATTCAC/p.P643_P644insHS	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548159_23548160insGTGAAT;unalign_gDNA=g.23548159_23548160
   insGTGAAT;insertion_gDNA=AATGTG;left_align_cDNA=c.1927_1928insCACATT;unalign_
   cDNA=c.1930_1931insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.P643_P64
   4insHS;unalign_protein=p.P643_P644insHS;phase=1
ACIN1:c.1930_1931insATTCAC	CCDS55905 (protein_coding)	ACIN1	-
   chr14:g.23548792_23548793insTGTGAA/c.1930_1931insATTCAC/p.S643_R644insHS	inside_[cds_in_exon_6]
   left_align_gDNA=g.23548787_23548788insGTGAAT;unalign_gDNA=g.23548787_23548788
   insGTGAAT;insertion_gDNA=TGTGAA;left_align_cDNA=c.1925_1926insTTCACA;unalign_
   cDNA=c.1930_1931insATTCAC;insertion_cDNA=ATTCAC;left_align_protein=p.R642_S64
   3insSH;unalign_protein=p.S643_R644insHS;phase=1
```
Reverse annotation can result in different identifiers after left/right alignments, e.g., 
```
#!bash
$ transvar canno --ccds -i 'AATK:c.3976_3977insCGCCCA'
```
results in
```
AATK:c.3976_3977insCGCCCA	CCDS45807 (protein_coding)	AATK	-
   chr17:g.79093282_79093287dupTGGGCG/c.3988_3993dupACGCCC/p.T1330_P1331dupTP	inside_[cds_in_exon_13]
   left_align_gDNA=g.79093270_79093271insGGGCGT;unalign_gDNA=g.79093282_79093287
   dupTGGGCG;insertion_gDNA=TGGGCG;left_align_cDNA=c.3976_3977insCGCCCA;unalign_
   cDNA=c.3976_3977insCGCCCA;insertion_cDNA=ACGCCC;left_align_protein=p.A1326_P1
   327insPT;unalign_protein=p.A1326_P1327insPT;phase=1
```
Note how insertion switch to duplication when 5'flanking is identical. This conforms to HGVS recommendation to replace insertion notation with duplication when possible.

Example: to annotate a **frame-shift insertion**, frameshift mutations have not alternative alignments. Hence only cDNA and gDNA have left alignment and unalignment reports.
```
#!bash
$ transvar canno --ccds -i 'AAAS:c.1225_1226insG'
```
results in
```
#!text
AAAS:c.1225_1226insG	CCDS8856 (protein_coding)	AAAS	-
   chr12:g.53702093dupC/c.1225dupG/p.E409Gfs*17	inside_[cds_in_exon_13]
   left_align_gDNA=g.53702089_53702090insC;unalign_gDNA=g.53702089_53702090insC;
   insertion_gDNA=C;left_align_cDNA=c.1221_1222insG;unalign_cDNA=c.1225dupG;inse
   rtion_cDNA=G
AAAS:c.1225_1226insG	CCDS53797 (protein_coding)	AAAS	-
   chr12:g.53701842_53701843insC/c.1225_1226insG/p.L409Rfs*54	inside_[cds_in_exon_13]
   left_align_gDNA=g.53701842_53701843insC;unalign_gDNA=g.53701842_53701843insC;
   insertion_gDNA=C;left_align_cDNA=c.1225_1226insG;unalign_cDNA=c.1225_1226insG
   ;insertion_cDNA=G
```

Example: to annotate an **intronic insertion**,
```
#!bash
$ transvar canno --ccds -i 'ADAM33:c.991-3_991-2insC'
```
outputs
```
#!text
ADAM33:c.991-3_991-2insC	CCDS13058 (protein_coding)	ADAM33	-
   chr20:g.3654151dupG/c.991-3dupC/.	inside_[intron_between_exon_10_and_11]
   left_align_gDNA=g.3654145_3654146insG;unalign_gDNA=g.3654145_3654146insG;inse
   rtion_gDNA=G;left_align_cDNA=c.991-9_991-8insC;unalign_cDNA=c.991-3dupC;inser
   tion_cDNA=C
```
In the case of intronic insertions, amino acid identifier is not applicable, represented in a `.`. But cDNA and gDNA identifier are right-aligned according to their natural order, respecting HGVS nomenclature.

Insertion could occur to *splice sites*. TransVar identifies such cases and report splice site and repress translation of protein change.
```
#!bash
$ transvar canno --ccds -i 'ADAM33:c.991_992insC'
```
results in
```
#!text
ADAM33:c.991_992insC	CCDS13058 (protein_coding)	ADAM33	-
   chr20:g.3654142_3654143insG/c.991_992insC/.	inside_[cds_in_exon_11]
   left_align_gDNA=g.3654142_3654143insG;unalign_gDNA=g.3654142_3654143insG;inse
   rtion_gDNA=G;left_align_cDNA=c.991_992insC;unalign_cDNA=c.991_992insC;inserti
   on_cDNA=C;acceptor_splice_site_on_exon_11_at_chr20:3654144
```

---

#### reverse annotation of nucleotide deletion
Similar to insertions, deletion can be in-frame or frame-shift. The consequence of deletion to amino acid sequence may appear a simple deletion or a block substitution (in the case where in-frame deletion is out of phase, i.e., partially delete codons).

Example: to annotate an **in-frame deletion**,
```
#!bash
$ transvar canno --ccds -i 'A4GNT:c.694_696delTTG'
```
```
#!text
A4GNT:c.694_696delTTG	CCDS3097 (protein_coding)	A4GNT	-
   chr3:g.137843435_137843437delACA/c.694_696delTTG/p.L232delL	inside_[cds_in_exon_2]
   left_align_gDNA=g.137843433_137843435delCAA;unaligned_gDNA=g.137843433_137843
   435delCAA;left_align_cDNA=c.692_694delTGT;unalign_cDNA=c.694_696delTTG;left_a
   lign_protein=p.L232delL;unalign_protein=p.L232delL;deletion_gDNA=ACA;deletion
   _cDNA=TTG
```

Example: to annotate a **in-frame, out-of-phase deletion**,
```
#!bash
$ transvar canno --ccds -i 'ABHD15:c.431_433delGTG'
```
```
#!text
ABHD15:c.431_433delGTG	CCDS32602 (protein_coding)	ABHD15	-
   chr17:g.27893552_27893554delCAC/c.431_433delGTG/p.C144_V145delinsF	inside_[cds_in_exon_1]
   left_align_gDNA=g.27893552_27893554delCAC;unaligned_gDNA=g.27893552_27893554d
   elCAC;left_align_cDNA=c.431_433delGTG;unalign_cDNA=c.431_433delGTG;deletion_g
   DNA=CAC;deletion_cDNA=GTG
```

Example: to annotate a **frame-shift deletion**,
```
#!bash
$ transvar canno --ccds -i 'AADACL3:c.374delG'
```
```
#!text
AADACL3:c.374delG	CCDS41252 (protein_coding)	AADACL3	+
   chr1:g.12785494delG/c.374delG/p.C125Ffs*17	cds_in_exon_3
   left_align_gDNA=g.12785494delG;unaligned_gDNA=g.12785494delG;left_align_cDNA=
   c.374delG;unalign_cDNA=c.374delG;deletion_gDNA=G;deletion_cDNA=G
```

Example: to annotate a **deletion that span from intronic to coding region**, protein prediction is suppressed due to loss of splice site.
```
#!bash
$ transvar canno --ccds -i 'ABCB11:c.1198-8_1199delcactccagAA'
```
```
#!text
ABCB11:c.1198-8_1199delcactccagAA	CCDS46444 (protein_coding)	ABCB11	-
   chr2:g.169833196_169833205delTTCTGGAGTG/c.1198-8_1199delCACTCCAGAA/.	from_[cds_in_exon_11]_to_[intron_between_exon_10_and_11]
   left_align_gDNA=g.169833196_169833205delTTCTGGAGTG;unaligned_gDNA=g.169833196
   _169833205delTTCTGGAGTG;left_align_cDNA=c.1198-8_1199delCACTCCAGAA;unalign_cD
   NA=c.1198-8_1199delCACTCCAGAA;acceptor_splice_site_on_exon_11_at_chr2:1698331
   98_lost;deletion_gDNA=TTCTGGAGTG;deletion_cDNA=CACTCCAGAA
```

---

#### reverse annotation of nucleotide block substitution

Example: to annotate a block substitution in **coding region**,
```
#!bash
$ transvar canno --ccds -i 'A1CF:c.508_509delinsTT'
```
```
#!text
A1CF:c.508_509delinsTT	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
A1CF:c.508_509delinsTT	CCDS7242 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
A1CF:c.508_509delinsTT	CCDS7243 (protein_coding)	A1CF	-
   chr10:g.52595953_52595954delinsAA/c.508_509delinsTT/p.G170F	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
```

Block substitution does not necessarily results in block substitution in amino acid. For example, the following substitution results in a deletion, where protein alternative alignment should be reported.
```
#!bash
$ transvar canno --ccds -i 'CSRNP1:c.1212_1224delinsGGAGGAGGAA'
```
```
#!text
CSRNP1:c.1212_1224delinsGGAGGAGGAA	CCDS2682 (protein_coding)	CSRNP1	-
   chr3:g.39185092_39185104delinsTTCCTCCTCC/c.1212_1224delinsGGAGGAGGAA/p.E411delE	inside_[cds_in_exon_4]
   begin_codon_cDNA=1210-1211-1212;end_codon_cDNA=1222-1223-1224;left_align_prot
   ein=p.E405delE;unalign_protein=p.E408delE
```

Likewise, block substitution could occur to **intronic region**,
```
#!bash
$ transvar canno --ccds -i 'A1CF:c.1460+2_1460+3delinsCC'
```
```
#!text
A1CF:c.1460+2_1460+3delinsCC	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52570797_52570798delinsGG/c.1460+2_1460+3delinsCC/.	inside_[intron_between_exon_9_and_10]
   .
```

When block substitution occurs **across splice site**, TransVar put a tag in the info fields and does not predict amino acid change.
```
#!bash
$ transvar canno --ccds -i 'A1CF:c.1459_1460+3delinsCC'
```
```
#!text
A1CF:c.1459_1460+3delinsCC	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52570797_52570801delinsGG/c.1459_1460+3delinsCC/.	from_[intron_between_exon_9_and_10]_to_[cds_in_exon_9]
   donor_splice_site_on_exon_9_at_chr10:52570799
```

---

#### reverse annotation of nucleotide duplication

Duplication can be thought of as special insertion where the inserted sequence is identical to the sequence flanking the breakpoint.
Similar to insertion, the annotation of duplication may possess alternative alignment.

Example: to annotate a duplication coding region,
```
#!bash
$ transvar canno --ccds -i 'CHD7:c.1669_1674dup'
```
```
#!text
CHD7:c.1669_1674dup	CCDS47865 (protein_coding)	CHD7	+
   chr8:g.61693564_61693569dupCCCGTC/c.1669_1674dup/p.P558_S559dupPS	inside_[cds_in_exon_2]
   left_align_gDNA=g.61693561_61693562insTCCCCG;unalign_gDNA=g.61693562_61693567
   dupTCCCCG;insertion_gDNA=CCCGTC;left_align_cDNA=c.1668_1669insTCCCCG;unalign_
   cDNA=c.1669_1674dupTCCCCG;insertion_cDNA=CCCGTC;left_align_protein=p.H556_S55
   7insSP;unalign_protein=p.S557_P558dupSP;phase=0
```

Example: a duplication on the nucleotide level may lead to frame-shift or block substitution on the amino acid level,
```
#!bash
$ transvar canno --ccds -i 'CHD7:c.1668_1669dup'
```
```
#!text
CHD7:c.1668_1669dup	CCDS47865 (protein_coding)	CHD7	+
   chr8:g.61693561_61693562dupTT/c.1668_1669dup/p.S557Ffs*8	inside_[cds_in_exon_2]
   left_align_gDNA=g.61693560_61693561insTT;unalign_gDNA=g.61693561_61693562dupT
   T;insertion_gDNA=TT;left_align_cDNA=c.1667_1668insTT;unalign_cDNA=c.1668_1669
   dupTT;insertion_cDNA=TT
```

Example: to annotate a duplication in intronic region,
```
#!bash
$ transvar canno --ccds -i 'CHD7:c.1666-5_1666-3dup'
```
```
#!text
CHD7:c.1666-5_1666-3dup	CCDS47865 (protein_coding)	CHD7	+
   chr8:g.61693554_61693556dupCTC/c.1666-5_1666-3dup/.	inside_[intron_between_exon_1_and_2]
   left_align_gDNA=g.61693553_61693554insCTC;unalign_gDNA=g.61693554_61693556dup
   CTC;insertion_gDNA=CTC;left_align_cDNA=c.1666-6_1666-5insCTC;unalign_cDNA=c.1
   666-5_1666-3dupCTC;insertion_cDNA=CTC
```

---

#### reverse annotation of amino acid insertion

```
#!bash
$ transvar panno --ccds -i 'AATK:p.P1331_A1332insTP'
```
```
#!text
AATK:p.P1331_A1332insTP	CCDS45807 (protein_coding)	AATK	-
   chr17:g.(79093267ins6)/c.(3997_3991ins6)/p.T1330_P1331dupTP	cds_in_exon_13
   left_align_protein=p.A1326_P1327insPT;unalign_protein=p.T1330_P1331dupTP;inse
   rtion_cDNA=ACACCT;insertion_gDNA=AGGTGT;imprecise
```

#### reverse annotation of amino acid deletion
```
#!bash
$ transvar panno --ccds -i 'AADACL4:p.W263_I267delWRDAI'
```
```
#!text
AADACL4:p.W263_I267delWRDAI	CCDS30590 (protein_coding)	AADACL4	+
   chr1:g.12726309_12726323del/c.787_801del/p.W263_I267delWRDAI	inside_[cds_in_exon_4]
   left_align_protein=p.W263_I267delWRDAI;unalign_protein=p.W263_I267delWRDAI;im
   precise
```

#### reverse annotation of amino acid block substitution
```
#!bash
$ transvar panno --ccds -i 'ABCC3:p.Y556_V557delinsRRR'
```
```
#!text
ABCC3:p.Y556_V557delinsRRR	CCDS32681 (protein_coding)	ABCC3	+
   chr17:g.48745254_48745259delinsAGGAGGAGG/c.1666_1671delinsAGGAGGAGG/p.Y556_V557delinsRRR	cds_in_exon_13
   imprecise
```

#### reverse annotation of amino acid frame-shift

```
#!bash
$ transvar panno --ccds -i 'A1BG:p.G132fs*2'
```
```
#!text
A1BG:p.G132fs*2	CCDS12976 (protein_coding)	A1BG	-
   chr19:g.58863860-58863868/c.394-402/p.G132fs*2	cds_in_exon_4
   imprecise
A1BG:p.G132fs*2	CCDS12976 (protein_coding)	A1BG	-
   chr19:g.58863860-58863868/c.394-402/p.G132fs*2	cds_in_exon_4
   imprecise
```

---

#### search alternative codon identifiers

An identifier is regarded as an alternative if the underlying codon overlap with the one from the original identifier.
Example: to search alternative identifiers of CDKN2A.p.58 (without knowing reference allele),
```
#!bash
$ transvar codonsearch --ccds -i CDKN2A:p.58
```
```
#!text
origin_id	alt_id	chrm	codon1
   codon2	transcripts_choice
CDKN2A:p.58	CDKN2A.p.73	chr9	21971184-21971185-21971186
   21971182-21971183-21971184	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
CDKN2A:p.58	CDKN2A.p.72	chr9	21971184-21971185-21971186
   21971185-21971186-21971187	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
```
The pair of transcript id listed corresponds to the transcripts based on which, the original and alternative identifiers are defined. Multiple pairs of transcript definitions are appended following a `,`.

Example: to search alternative identifiers of DHODH:G152R (knowing reference allele `G`, alternative allele here will be ignored),
```
#!bash
$ transvar codonsearch -i DHODH:G152R --refseq
```
outputs
```
#!text
origin_id	alt_id	chrm	codon1
   codon2	transcripts_choice
DHODH:G152R	DHODH.p.G16	chr16	72050942-72050943-72050944
   72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255828[RefSeq]
DHODH:G152R	DHODH.p.G9	chr16	72050942-72050943-72050944
   72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255829[RefSeq]
DHODH:G152R	DHODH.p.G124	chr16	72050942-72050943-72050944
   72050942-72050943-72050944	NM_001361[RefSeq]/XM_005255827[RefSeq]
```
TransVar outputs genomic positions of codons based on original transcript (4th column in the output) and alternative transcript (5th column in the output). The potential transcript usages are also appended.

Example: to run `transvar codonsearch` to **batch process** a list of mutation identifiers.
```
#!bash
$ transvar codonsearch -l example/input_table2 --ccds -m 1 -o 1
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
origin_id	alt_id	chrm	codon1
   codon2	transcripts_choice
CDKN2A:p.61	CDKN2A.p.76	chr9	21971175-21971176-21971177
   21971173-21971174-21971175	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
CDKN2A:p.61	CDKN2A.p.75	chr9	21971175-21971176-21971177
   21971176-21971177-21971178	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
CDKN2A:p.69	CDKN2A.p.54	chr9	21971194-21971195-21971196
   21971196-21971197-21971198	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
CDKN2A:p.69	CDKN2A.p.55	chr9	21971194-21971195-21971196
   21971193-21971194-21971195	CCDS6511[CCDS]/CCDS6510[CCDS],CCDS6511[CCDS]/CCDS56565[CCDS]
CDKN2A:p.69	CDKN2A.p.83	chr9	21971151-21971152-21971153
   21971152-21971153-21971154	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
CDKN2A:p.69	CDKN2A.p.84	chr9	21971151-21971152-21971153
   21971149-21971150-21971151	CCDS6510[CCDS]/CCDS6511[CCDS],CCDS56565[CCDS]/CCDS6511[CCDS]
ERBB2:p.755	ERBB2.p.785	chr17	37881024-37881025-37881026
   37881024-37881025-37881026	CCDS45667[CCDS]/CCDS32642[CCDS]
ERBB2:p.755	ERBB2.p.725	chr17	37880219-37880220-37880221
   37880219-37880220-37880221	CCDS32642[CCDS]/CCDS45667[CCDS]
```
The third column indicates the potential transcript usage for the alternative identifier. Each transcript usage is denoted by <listing transcript>/<actual transcript>. Different potential choices are separated by ','.

---
#### infer potential codon identity

Example: to check if MET.p1010 and MET.p992 may be refering to one mutation due to different usage of transcripts,
```
#!bash
$ transvar codonsearch --refseq -i MET:p.1010
```
gives
```
#!text
origin_id	alt_id	chrm	codon1
   codon2	transcripts_choice
MET:p.1010	MET.p.562	chr7	116411989-116411990-116411991
   116411989-116411990-116411991	NM_001127500[RefSeq]/XM_005250354[RefSeq]
MET:p.1010	MET.p.1029	chr7	116411989-116411990-116411991
   116411989-116411990-116411991	NM_001127500[RefSeq]/XM_005250353[RefSeq]
MET:p.1010	MET.p.973	chr7	116411932-116411933-116411934
   116411932-116411933-116411934	XM_005250353[RefSeq]/NM_000245[RefSeq]
MET:p.1010	MET.p.580	chr7	116412043-116414935-116414936
   116412043-116414935-116414936	NM_000245[RefSeq]/XM_005250354[RefSeq]
MET:p.1010	MET.p.991	chr7	116411932-116411933-116411934
   116411932-116411933-116411934	XM_005250353[RefSeq]/NM_001127500[RefSeq]
MET:p.1010	MET.p.543	chr7	116411932-116411933-116411934
   116411932-116411933-116411934	XM_005250353[RefSeq]/XM_005250354[RefSeq]
MET:p.1010	MET.p.1028	chr7	116412043-116414935-116414936
   116412043-116414935-116414936	NM_000245[RefSeq]/NM_001127500[RefSeq]
MET:p.1010	MET.p.992	chr7	116411989-116411990-116411991
   116411989-116411990-116411991	NM_001127500[RefSeq]/NM_000245[RefSeq]
MET:p.1010	MET.p.1047	chr7	116412043-116414935-116414936
   116412043-116414935-116414936	NM_000245[RefSeq]/XM_005250353[RefSeq]
```
Since MET.p.992 is in the list, the two identifiers might be due to the same genomic mutation.

#### annotate SNP from genomic locations

This is the forward annotation

```
#!bash
$ transvar ganno --ccds -i 'chr3:g.178936091G>A'
```
outputs
```
#!text
chr3:g.178936091G>A	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	cds_in_exon_9
   dbsnp=rs104886003(chr3:178936091G>A);missense;codon_pos=178936091-178936092-1
   78936093;ref_codon_seq=GAG
```

Another example:
```
#!bash
$ transvar ganno -i "chr9:g.135782704C>G" --ccds
```
outputs
```
#!text
chr9:g.135782704C>G	CCDS6956 (protein_coding)	TSC1	-
   chr9:g.135782704C>G/c.1317G>C/p.L439L	cds_in_exon_11
   dbsnp=rs770692313(chr9:135782704C>G);synonymous;codon_pos=135782704-135782705
   -135782706;ref_codon_seq=CTG
chr9:g.135782704C>G	CCDS55350 (protein_coding)	TSC1	-
   chr9:g.135782704C>G/c.1164G>C/p.L388L	cds_in_exon_10
   dbsnp=rs770692313(chr9:135782704C>G);synonymous;codon_pos=135782704-135782705
   -135782706;ref_codon_seq=CTG
```


#### annotate a short genomic region

To annotate a short genomic region in a gene,
```
#!bash
$ transvar ganno --ccds -i 'chr3:g.178936091_178936192'
```
outputs
```
#!text
chr3:g.178936091_178936192	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091_178936192/c.1633_1664+70/p.E545_R555	from_[cds_in_exon_9]_to_[intron_between_exon_9_and_10]
   donor_splice_site_on_exon_9_at_chr3:178936123_included;start_codon=178936091-
   178936092-178936093;end_codon=178936121-178936122-178936984
```
	
Results indicates the beginning position is at coding region while ending position is at intronic region (c.1633_1664+70).

For intergenic sites, TransVar also reports the identity and distance to the gene upstream and downstream. For example, `chr6:116991832` is simply annotated as intergenic in the original annotation. TransVar reveals that it is 1,875 bp downstream to ZUFSP and 10,518 bp upstream to KPNA5 showing a vicinity to the gene ZUFSP. There is no limit in the reported distance. If a site is at the end of the chromosome, TransVar is able to report the distance to the telomere.

#### annotate a long genomic region
[back to top](#top)
```
#!bash
$ transvar ganno -i '9:g.133750356_137990357' --ccds
```
outputs
```
#!text
9:g.133750356_137990357	CCDS35165 (protein_coding),CCDS6986 (protein_coding)	.	.
   chr9:g.133750356_137990357/./.	from_[cds_in_exon_7;ABL1]_to_[intron_between_exon_4_and_5;OLFM1]_spanning_[51_genes]
   .
9:g.133750356_137990357	CCDS35166 (protein_coding),CCDS6986 (protein_coding)	.	.
   chr9:g.133750356_137990357/./.	from_[cds_in_exon_7;ABL1]_to_[intron_between_exon_4_and_5;OLFM1]_spanning_[51_genes]
   .
```
The result indicates that the region span 53 genes. The beginning of the region resides in the coding sequence of ABL1, c.1187A and the ending region resides in the intronic region of OLFM1, c.622+6C. 2 different usage of transcripts in annotating the starting position is represented in two lines, each line corresponding to a combination of transcript usage.
This annotation not only shows the coverage of the region, also reveals the fine structure of the boundary.

In another example, where the ending position exceeds the length of the chromosome, TransVar truncates the region and outputs upstream and downstream information of the ending position.
```
#!bash
$ transvar ganno -i '9:g.133750356_1337503570' --ccds
```
outputs
```
#!text
9:g.133750356_1337503570	CCDS35165 (protein_coding),	.	.
   chr9:g.133750356_141213431/./.	from_[cds_in_exon_7;ABL1]_to_[intergenic_between_EHMT1(484,026_bp_downstream)_and_3'-telomere(0_bp)]_spanning_[136_genes]
   .
9:g.133750356_1337503570	CCDS35166 (protein_coding),	.	.
   chr9:g.133750356_141213431/./.	from_[cds_in_exon_7;ABL1]_to_[intergenic_between_EHMT1(484,026_bp_downstream)_and_3'-telomere(0_bp)]_spanning_[136_genes]
   .
```

#### annotate a deletion from genomic location
[back to top](#top)

A frameshift deletion
```
#!bash
$ transvar ganno -i "chr2:g.234183368_234183380del" --ccds
```
outputs
```
#!text
chr2:g.234183368_234183380del	CCDS2502 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183380del13/c.841_853del13/p.T281Lfs*5	inside_[cds_in_exon_8]
   left_align_gDNA=g.234183367_234183379del13;unaligned_gDNA=g.234183368_2341833
   80del13;left_align_cDNA=c.840_852del13;unalign_cDNA=c.841_853del13
chr2:g.234183368_234183380del	CCDS2503 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183380del13/c.898_910del13/p.T300Lfs*5	inside_[cds_in_exon_9]
   left_align_gDNA=g.234183367_234183379del13;unaligned_gDNA=g.234183368_2341833
   80del13;left_align_cDNA=c.897_909del13;unalign_cDNA=c.898_910del13
chr2:g.234183368_234183380del	CCDS54438 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183380del13/c.409_421del13/p.T137Lfs*5	inside_[cds_in_exon_5]
   left_align_gDNA=g.234183367_234183379del13;unaligned_gDNA=g.234183368_2341833
   80del13;left_align_cDNA=c.408_420del13;unalign_cDNA=c.409_421del13
```
Note the difference between left-aligned identifier and the right aligned identifier.

An in-frame deletion
```
#!bash
$ transvar ganno -i "chr2:g.234183368_234183379del" --ccds
```
outputs
```
#!text
chr2:g.234183368_234183379del	CCDS2502 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183379del12/c.841_852del12/p.T281_G284delTHPG	inside_[cds_in_exon_8]
   left_align_gDNA=g.234183367_234183378del12;unaligned_gDNA=g.234183368_2341833
   79del12;left_align_cDNA=c.840_851del12;unalign_cDNA=c.841_852del12;left_align
   _protein=p.T281_G284delTHPG;unalign_protein=p.T281_G284delTHPG
chr2:g.234183368_234183379del	CCDS2503 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183379del12/c.898_909del12/p.T300_G303delTHPG	inside_[cds_in_exon_9]
   left_align_gDNA=g.234183367_234183378del12;unaligned_gDNA=g.234183368_2341833
   79del12;left_align_cDNA=c.897_908del12;unalign_cDNA=c.898_909del12;left_align
   _protein=p.T300_G303delTHPG;unalign_protein=p.T300_G303delTHPG
chr2:g.234183368_234183379del	CCDS54438 (protein_coding)	ATG16L1	+
   chr2:g.234183368_234183379del12/c.409_420del12/p.T137_G140delTHPG	inside_[cds_in_exon_5]
   left_align_gDNA=g.234183367_234183378del12;unaligned_gDNA=g.234183368_2341833
   79del12;left_align_cDNA=c.408_419del12;unalign_cDNA=c.409_420del12;left_align
   _protein=p.T137_G140delTHPG;unalign_protein=p.T137_G140delTHPG
```

Another example
```
#!bash
$ transvar ganno --ccds -i 'chr12:g.53703425_53703427del'
```
outputs
```
#!text
chr12:g.53703425_53703427del	CCDS8856 (protein_coding)	AAAS	-
   chr12:g.53703427_53703429delCCC/c.769_771delGGG/p.G257delG	inside_[cds_in_exon_8]
   left_align_gDNA=g.53703424_53703426delCCC;unaligned_gDNA=g.53703425_53703427d
   elCCC;left_align_cDNA=c.766_768delGGG;unalign_cDNA=c.768_770delGGG;left_align
   _protein=p.G256delG;unalign_protein=p.G256delG
chr12:g.53703425_53703427del	CCDS53797 (protein_coding)	AAAS	-
   chr12:g.53703427_53703429delCCC/c.670_672delGGG/p.G224delG	inside_[cds_in_exon_7]
   left_align_gDNA=g.53703424_53703426delCCC;unaligned_gDNA=g.53703425_53703427d
   elCCC;left_align_cDNA=c.667_669delGGG;unalign_cDNA=c.669_671delGGG;left_align
   _protein=p.G223delG;unalign_protein=p.G223delG
```
Note the difference between left and right-aligned identifiers on both protein level and cDNA level.

An in-frame out-of-phase deletion
```
#!bash
$ transvar ganno -i "chr2:g.234183372_234183383del" --ccds
```
outputs
```
#!text
chr2:g.234183372_234183383del	CCDS2502 (protein_coding)	ATG16L1	+
   chr2:g.234183372_234183383del12/c.845_856del12/p.H282_G286delinsR	inside_[cds_in_exon_8]
   left_align_gDNA=g.234183372_234183383del12;unaligned_gDNA=g.234183372_2341833
   83del12;left_align_cDNA=c.845_856del12;unalign_cDNA=c.845_856del12
chr2:g.234183372_234183383del	CCDS2503 (protein_coding)	ATG16L1	+
   chr2:g.234183372_234183383del12/c.902_913del12/p.H301_G305delinsR	inside_[cds_in_exon_9]
   left_align_gDNA=g.234183372_234183383del12;unaligned_gDNA=g.234183372_2341833
   83del12;left_align_cDNA=c.902_913del12;unalign_cDNA=c.902_913del12
chr2:g.234183372_234183383del	CCDS54438 (protein_coding)	ATG16L1	+
   chr2:g.234183372_234183383del12/c.413_424del12/p.H138_G142delinsR	inside_[cds_in_exon_5]
   left_align_gDNA=g.234183372_234183383del12;unaligned_gDNA=g.234183372_2341833
   83del12;left_align_cDNA=c.413_424del12;unalign_cDNA=c.413_424del12
```

#### annotate an insertion from genomic location

An in-frame insertion of three nucleotides
```
#!bash
$ transvar ganno -i 'chr2:g.69741762_69741763insTGC' --ccds
```
outputs
```
#!text
chr2:g.69741762_69741763insTGC	CCDS1893 (protein_coding)	AAK1	-
   chr2:g.69741780_69741782dupCTG/c.1614_1616dupGCA/p.Q546dupQ	cds_in_exon_12
   left_align_gDNA=g.69741762_69741763insTGC;unalign_gDNA=g.69741762_69741763ins
   TGC;insertion_gDNA=CTG;left_align_cDNA=c.1596_1597insCAG;unalign_cDNA=c.1614_
   1616dupGCA;insertion_cDNA=GCA;left_align_protein=p.Y532_Q533insQ;unalign_prot
   ein=p.Q539dupQ;phase=2
```
Note the proper right-alignment of protein level insertion Q. The left-aligned identifier is also given in the `LEFTALN` field.

A frame-shift insertion of two nucleotides
```
#!bash
$ transvar ganno -i 'chr7:g.121753754_121753755insCA' --ccds
```
outputs
```
#!text
chr7:g.121753754_121753755insCA	CCDS5783 (protein_coding)	AASS	-
   chr7:g.121753754_121753755insCA/c.1064_1065insGT/p.I355Mfs*10	cds_in_exon_9
   left_align_gDNA=g.121753753_121753754insAC;unalign_gDNA=g.121753754_121753755
   insCA;insertion_gDNA=CA;left_align_cDNA=c.1063_1064insTG;unalign_cDNA=c.1063_
   1064insTG;insertion_cDNA=GT
```

```
#!bash
$ transvar ganno -i 'chr17:g.79093270_79093271insGGGCGT' --ccds
```
outputs
```
#!text
chr17:g.79093270_79093271insGGGCGT	CCDS45807 (protein_coding)	AATK	-
   chr17:g.79093282_79093287dupTGGGCG/c.3988_3993dupACGCCC/p.T1330_P1331dupTP	cds_in_exon_13
   left_align_gDNA=g.79093270_79093271insGGGCGT;unalign_gDNA=g.79093270_79093271
   insGGGCGT;insertion_gDNA=TGGGCG;left_align_cDNA=c.3976_3977insCGCCCA;unalign_
   cDNA=c.3988_3993dupACGCCC;insertion_cDNA=ACGCCC;left_align_protein=p.A1326_P1
   327insPT;unalign_protein=p.T1330_P1331dupTP;phase=0
```
Notice the difference in the inserted sequence when left-alignment and right-alignment conventions are followed.

A frame-shift insertion of one nucleotides in a homopolymer
```
#!bash
$ transvar ganno -i 'chr7:g.117230474_117230475insA' --ccds
```
outputs
```
#!text
chr7:g.117230474_117230475insA	CCDS5773 (protein_coding)	CFTR	+
   chr7:g.117230479dupA/c.1752dupA/p.E585Rfs*4	cds_in_exon_13
   left_align_gDNA=g.117230474_117230475insA;unalign_gDNA=g.117230474_117230475i
   nsA;insertion_gDNA=A;left_align_cDNA=c.1747_1748insA;unalign_cDNA=c.1747_1748
   insA;insertion_cDNA=A
```
Notice the right alignment of cDNA level insertion and the left alignment reported as additional information.

A in-frame, in-phase insertion
```
#!bash
$ transvar ganno -i 'chr12:g.109702119_109702120insACC' --ccds
```
```
#!text
chr12:g.109702119_109702120insACC	CCDS31898 (protein_coding)	ACACB	+
   chr12:g.109702119_109702120insACC/c.6870_6871insACC/p.Y2290_H2291insT	cds_in_exon_49
   left_align_gDNA=g.109702118_109702119insCAC;unalign_gDNA=g.109702119_10970212
   0insACC;insertion_gDNA=ACC;left_align_cDNA=c.6869_6870insCAC;unalign_cDNA=c.6
   870_6871insACC;insertion_cDNA=ACC;left_align_protein=p.Y2290_H2291insT;unalig
   n_protein=p.Y2290_H2291insT;phase=0
```

#### annotate block substitution from genomic locations

A block-substitution that results in a frameshift.
```
#!bash
$ transvar ganno -i 'chr10:g.27329002_27329002delinsAT' --ccds
```
```
#!text
chr10:g.27329002_27329002delinsAT	CCDS41499 (protein_coding)	ANKRD26	-
   chr10:g.27329009dupT/c.2266dupA/p.M756Nfs*6	cds_in_exon_21
   left_align_gDNA=g.27329002_27329003insT;unalign_gDNA=g.27329002_27329003insT;
   insertion_gDNA=T;left_align_cDNA=c.2259_2260insA;unalign_cDNA=c.2266dupA;inse
   rtion_cDNA=A
```

A block-substitution that is in-frame,
```
#!bash
$ transvar ganno -i 'chr10:g.52595929_52595930delinsAA' --ccds
```
```
#!text
chr10:g.52595929_52595930delinsAA	CCDS7243 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.532_533delinsTT/p.P178L	inside_[cds_in_exon_4]
   codon_cDNA=532-533-534
chr10:g.52595929_52595930delinsAA	CCDS7242 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
chr10:g.52595929_52595930delinsAA	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
```

#### annotate promoter region

One can define the promoter boundary through the `--prombeg` and `--promend` option. Default promoter region is defined from 1000bp upstream of the transcription start site to the transcription start site. One could customize this setting to e.g., [-1000bp, 2000bp] by

```
#!bash
$ transvar ganno -i 'chr19:g.41950335_41951908' --ensembl --prombeg 2000 --promend 1000 --refversion mm10
```
```
#!text
chr19:g.41950335_41951908	ENSMUST00000167927 (nonsense_mediated_decay)	MMS19	-
   chr19:g.41950335_41951908/c.1071+3684_1071+5257/.	from_[intron_between_exon_20_and_21]_to_[intron_between_exon_19_and_20]
   whole_exon_[20]_included
chr19:g.41950335_41951908	ENSMUST00000171561 (protein_coding)	MMS19	-
   chr19:g.41950335_41951908/c.1915+499_2016-252/p.E639_E672	from_[intron_between_exon_20_and_21]_to_[intron_between_exon_19_and_20]
   whole_exon_[20]_included;start_codon=41950753-41950752-41950083;end_codon=419
   52407-41950851-41950850
chr19:g.41950335_41951908	ENSMUST00000170209 (retained_intron)	MMS19	-
   chr19:g.41950335_41951908/c.2251+499_2352-252/.	from_[intron_between_exon_16_and_17]_to_[intron_between_exon_15_and_16]
   whole_exon_[16]_included
chr19:g.41950335_41951908	ENSMUST00000163287 (protein_coding)	MMS19	-
   chr19:g.41950335_41951908/c.1477+499_1578-252/p.E493_E526	from_[intron_between_exon_17_and_18]_to_[intron_between_exon_16_and_17]
   whole_exon_[17]_included;start_codon=41950753-41950752-41950083;end_codon=419
   52407-41950851-41950850
chr19:g.41950335_41951908	ENSMUST00000163398 (nonsense_mediated_decay)	MMS19	-
   chr19:g.41950335_41951908/c.225+12487_225+14060/.	from_[intron_between_exon_19_and_20]_to_[intron_between_exon_18_and_19]
   whole_exon_[19]_included
chr19:g.41950335_41951908	ENSMUST00000164776 (nonsense_mediated_decay)	MMS19	-
   chr19:g.41950335_41951908/c.225+12487_225+14060/.	from_[intron_between_exon_19_and_20]_to_[intron_between_exon_18_and_19]
   whole_exon_[19]_included
chr19:g.41950335_41951908	ENSMUST00000026168 (protein_coding)	MMS19	-
   chr19:g.41950335_41951908/c.1786+499_1887-252/p.E596_E629	from_[intron_between_exon_19_and_20]_to_[intron_between_exon_18_and_19]
   whole_exon_[19]_included;start_codon=41950753-41950752-41950083;end_codon=419
   52407-41950851-41950850
chr19:g.41950335_41951908	ENSMUST00000166090 (nonsense_mediated_decay)	MMS19	-
   chr19:g.41950335_41951908/c.636+499_737-252/.	from_[intron_between_exon_7_and_8]_to_[intron_between_exon_6_and_7]
   whole_exon_[7]_included
chr19:g.41950335_41951908	ENSMUST00000171755 (retained_intron)	MMS19	-
   chr19:g.41950335_41951908/c.1941+499_2042-252/.	from_[intron_between_exon_20_and_21]_to_[intron_between_exon_19_and_20]
   whole_exon_[20]_included
chr19:g.41950335_41951908	ENSMUST00000167820 (protein_coding)	MMS19	-
   chr19:g.41950335_41951908/c.179-1057_279-252/p.E60_E93	from_[intron_between_exon_3_and_4]_to_[intron_between_exon_2_and_3]
   whole_exon_[3]_included;start_codon=41950753-41950752-41950083;end_codon=4195
   3669-41950851-41950850
chr19:g.41950335_41951908	ENSMUST00000169775 (nonsense_mediated_decay)	MMS19	-
   chr19:g.41950335_41951908/c.522+11101_522+12674/.	from_[intron_between_exon_20_and_21]_to_[intron_between_exon_19_and_20]
   whole_exon_[20]_included
chr19:g.41950335_41951908	ENSMUST00000166517 (retained_intron)	MMS19	-
   chr19:g.41950335_41951908/c.1-564_594-252/.	from_[intron_between_exon_1_and_2]_to_[intergenic_between_MMS19(564_bp_upstream)_and_MMS19(1,189_bp_downstream)]
   promoter_region_of_[MMS19]_overlaping_1565_bp(99.43%);whole_exon_[1]_included
```
The result shows that 99.43% of the target region is inside the promoter region. The overlap is as long as 1564 base pairs.

#### annotate non-coding RNA
Given Ensembl, GENCODE or RefSeq database, one could annotate non-coding transcripts such as lncRNA.
E.g.,
```
#!bash
$ transvar ganno --gencode -i 'chr1:g.3985200_3985300' --refversion mm10
```
results in
```
#!text
chr1:g.3985200_3985300	ENSMUST00000194643 (lincRNA)	RP23-333I7.1	-
   chr1:g.3985200_3985300/c.121_221/.	inside_[noncoding_exon_2]
   .
chr1:g.3985200_3985300	ENSMUST00000192427 (lincRNA)	RP23-333I7.1	-
   chr1:g.3985200_3985300/c.685_785/.	inside_[noncoding_exon_1]
   .
```
or
```
#!bash
$ transvar ganno --refseq -i 'chr14:g.20568338_20569581' --refversion mm10
```
results in
```
#!text
chr14:g.20568338_20569581	NR_033571 (lncRNA)	1810062O18RIK	+
   chr14:g.20568338_20569581/c.260-1532_260-289/.	inside_[intron_between_exon_4_and_5]
   dbxref=GeneID:75602,MGI:MGI:1922852
chr14:g.20568338_20569581	XM_011245228 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.1357+667_1357+1910/.	inside_[intron_between_exon_6_and_7]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_011245226 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.1972+667_1972+1910/.	inside_[intron_between_exon_13_and_14]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	NM_030180 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2188+667_2188+1910/.	inside_[intron_between_exon_15_and_16]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_011245225 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2359+667_2359+1910/.	inside_[intron_between_exon_16_and_17]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_006519705 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2188+667_2188+1910/.	inside_[intron_between_exon_15_and_16]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_006519703 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2359+667_2359+1910/.	inside_[intron_between_exon_16_and_17]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_011245227 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2359+667_2359+1910/.	inside_[intron_between_exon_16_and_17]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_006519709 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2359+667_2359+1910/.	inside_[intron_between_exon_16_and_17]
   dbxref=GeneID:78787,MGI:MGI:1926037
chr14:g.20568338_20569581	XM_006519708 (protein_coding)	USP54	-
   chr14:g.20568338_20569581/c.2359+667_2359+1910/.	inside_[intron_between_exon_16_and_17]
   dbxref=GeneID:78787,MGI:MGI:1926037
```

or using Ensembl
```
#!bash
$ transvar ganno --ensembl -i 'chr1:g.29560_29570'
```
results in
```
#!text
chr1:g.29560_29570	ENST00000488147 (unprocessed_pseudogene)	WASH7P	-
   chr1:g.29560_29570/c.1_11/.	inside_[noncoding_exon_1]
   promoter_region_of_[WASH7P]_overlaping_1_bp(9.09%)
chr1:g.29560_29570	ENST00000538476 (unprocessed_pseudogene)	WASH7P	-
   chr1:g.29560_29570/c.237_247/.	inside_[noncoding_exon_1]
   .
chr1:g.29560_29570	ENST00000473358 (lincRNA)	MIR1302-10	+
   chr1:g.29560_29570/c.7_17/.	inside_[noncoding_exon_1]
   .
```

### FAQ

#### how can TransVar take VCF as input?

Yes. For example,
```
#!bash
transvar ganno --vcf ALL.wgs.phase1_release_v3.20101123.snps_indel_sv.sites.vcf.gz --ccds
# or
transvar ganno --vcf demo.1kg.vcf --ccds
```

#### Can TransVar automatically decompose a haplotype into multiple mutations?

Yes, TransVar performs local alignment to allow long haplotype to be decomposed into multiple mutations.

```
#!bash
$ transvar ganno --ccds -i 'chr20:g.645097_645111delinsGTGCGATACCCAGGAG' --haplotype
```
leads to 2 snv and one insertion
```
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645098G>T/c.141C>A/p.A47A	cds_in_exon_2
   synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350	cds_in_exon_2
   left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645102insA;insertio
   n_gDNA=A;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;insertion_c
   DNA=T
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	intron_between_exon_1_and_2
   acceptor_splice_site_of_exon_1_at_chr20:645106
```

#### Can TransVar use 3-letter code instead of 1-letter code for protein?

Yes, TransVar automatically infer whether the input is a 3-letter code or 1-letter code.
The output is default to 1-letter code. But can be switched to 3-letter code through the `--aa3` option.
For example,
```
#!bash
$ transvar panno --ccds -i 'PIK3CA:p.Glu545Lys' --aa3
```
```
PIK3CA:p.Glu545Lys	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.Glu545Lys	cds_in_exon_9
   reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_variants=chr3:g.17
   8936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G>A);missense
```

#### Can TransVar report results in one line for each query?

Yes, with `--oneline` option. This separates the outputs from each transcript by '|||'.

#### I got 'gene_not_recognized', what's wrong?

Most likely you forgot to specify a transcipt definition such as `--ccds` or `--ensembl`. Sometimes there are non-canonical names for genes, this can be fixed through the `--alias` option and specify an alias table. TransVar comes with alias table from UCSC knownGene.

#### Does TransVar support alternative format for MNV such as `c.508_509CC>TT`?

Yes, but only in input. For example, `c.508_509CC>TT`
```
#!bash
$ transvar canno --ccds -i 'A1CF:c.508_509CC>TT'
```
```
A1CF:c.508_509CC>TT	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
A1CF:c.508_509CC>TT	CCDS7242 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   codon_cDNA=508-509-510
```

#### Does TransVar support relaxed input without 'g.', 'c.' and 'p.'?

Yes, the 'g.', 'c.' and 'p.' are optional in the input. For example, `12:109702119insACC` is equally acceptable as `chr12:g.109702119_109702120insACC`. TransVar also accepts '>' in denoting MNV. E.g., `c.113G>TACTAGC` can be used in place of `c.113delGinsTACTAGC`. This is common in some database such as COSMIC.

#### When I annotate a variant for protein identifier, why would I end up getting results in another variant type?

TransVar follows in full the HGVS nomenclature while annotating protein level mutation identifiers. For example, a out-of-phase, in frame insertion, `ACIN1:c.1930_1931insATTCAC` will be annotated with `p.S643_R644insHS` rather than `R644delinsHSR`. Protein level mutation will be generated as if no nucleotide mutation information exists.

## Future work

 + add cytoband annotation
 + imprecise annotation
 + forward annotation of binding sites
 + forward annotation of structural variation breakpoints
 + begin codon and end codon in deletion
 + distinguish non-transcribable element and suppress promoter setting (like "retained intron")

## Bug report and feature request

Please direct any bugs to <zhouwanding@gmail.com>.

## Reference

submitted

## About
This work is a collaboration between Wanding Zhou, Tenghui Chen, Zechen Chong and Professor Ken Chen at UT MD Anderson Cancer Center.


