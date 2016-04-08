## How to batch-process?

For all mutation types, one can batch process a list of mutation identifiers with optional transcript id to constraint the search. Take SNV for example,
```bash
transvar panno -l example/input_table -g 1 -m 5 -t 2 --ensembl -o 2,3,4
```
As suggested by the command, TransVar takes as input the 1st column as gene and 4th column as identifier. The 2nd column will be used as the transcript id from Ensembl to constrain the alternative identifier search. The 2nd, 3rd and 5th columns are chosen to be output as a validation of TransVar's performance.

Input:
```text
ADAMTSL3        ENST00000286744 15:84442328     c.243G>A        p.W81*  Nonsense
ADAMTSL3        ENST00000286744 15:84442326     c.241T>C        p.W81R  Missense
ADAMTSL4        ENST00000369038 1:150530513     c.2270G>A       p.G757D Missense
ADCY2   ENST00000338316 5:7802364       c.2662G>A       p.V888I Missense
ADCY2   ENST00000338316 5:7802365       c.2663T>C       p.V888A Missense
```
Output:
```text
ENST00000286744|15:84442328|c.243G>A	ENST00000286744 (protein_coding)	ADAMTSL3	+
   chr15:g.84442327G>A/c.242G>A/p.W81*	cds_in_exon_4
   reference_codon=TGG;candidate_codons=TAA,TAG,TGA;candidate_snv_variants=chr15
   :g.84442328G>A;candidate_mnv_variants=chr15:g.84442327_84442328delGGinsAA;mis
   sense;aliases=ENSP00000286744;source=Ensembl
ENST00000286744|15:84442326|c.241T>C	ENST00000286744 (protein_coding)	ADAMTSL3	+
   chr15:g.84442326T>A/c.241T>A/p.W81R	cds_in_exon_4
   reference_codon=TGG;candidate_codons=AGG,AGA,CGA,CGC,CGG,CGT;candidate_snv_va
   riants=chr15:g.84442326T>C;candidate_mnv_variants=chr15:g.84442326_84442328de
   lTGGinsAGA,chr15:g.84442326_84442328delTGGinsCGA,chr15:g.84442326_84442328del
   TGGinsCGC,chr15:g.84442326_84442328delTGGinsCGT;missense;aliases=ENSP00000286
   744;source=Ensembl
ENST00000369038|1:150530513|c.2270G>A	ENST00000369038 (protein_coding)	ADAMTSL4	+
   chr1:g.150530513G>A/c.2270G>A/p.G757D	cds_in_exon_12
   reference_codon=GGT;candidate_codons=GAC,GAT;candidate_mnv_variants=chr1:g.15
   0530513_150530514delGTinsAC;missense;aliases=ENSP00000358034;source=Ensembl
ENST00000338316|5:7802364|c.2662G>A	ENST00000338316 (protein_coding)	ADCY2	+
   chr5:g.7802364G>A/c.2662G>A/p.V888I	cds_in_exon_21
   reference_codon=GTC;candidate_codons=ATC,ATA,ATT;candidate_mnv_variants=chr5:
   g.7802364_7802366delGTCinsATA,chr5:g.7802364_7802366delGTCinsATT;missense;ali
   ases=ENSP00000342952;source=Ensembl
ENST00000338316|5:7802365|c.2663T>C	ENST00000338316 (protein_coding)	ADCY2	+
   chr5:g.7802365T>C/c.2663T>C/p.V888A	cds_in_exon_21
   reference_codon=GTC;candidate_codons=GCA,GCC,GCG,GCT;candidate_mnv_variants=c
   hr5:g.7802365_7802366delTCinsCA,chr5:g.7802365_7802366delTCinsCG,chr5:g.78023
   65_7802366delTCinsCT;missense;aliases=ENSP00000342952;source=Ensembl
```

## How to use VCF as input?

TransVar can take VCF as input when annotating from genomic level.
```bash
transvar ganno --vcf ALL.wgs.phase1_release_v3.20101123.snps_indel_sv.sites.vcf.gz --ccds
# or
transvar ganno --vcf demo.1kg.vcf --ccds
```

## How to automatically decompose a haplotype into multiple mutations?

TransVar performs local alignment to allow long haplotype to be decomposed into multiple mutations.

```bash
$ transvar ganno --ccds -i 'chr20:g.645097_645111delinsGTGCGATACCCAGGAG' --haplotype
```
leads to 2 snv and one insertion
```text
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645098G>T/c.141C>A/p.A47A	inside_[cds_in_exon_2]
   CSQN=Synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC;source=CCDS
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350	inside_[cds_in_exon_2]
   CSQN=Frameshift;left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645
   102insA;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;source=CCDS
chr20:g.645097_645111delinsGTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	inside_[intron_between_exon_1_and_2]
   CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon1_At_chr20:645106;source=CCDS
```

## How to use 3-letter code instead of 1-letter code for protein?

TransVar automatically infer whether the input is a 3-letter code or 1-letter code.
The output is default to 1-letter code. But can be switched to 3-letter code through the `--aa3` option.
For example,
```bash
$ transvar panno --ccds -i 'PIK3CA:p.Glu545Lys' --aa3
```
```text
PIK3CA:p.Glu545Lys	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.Glu545Lys	inside_[cds_in_exon_9]
   CSQN=Missense;reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_vari
   ants=chr3:g.178936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G
   >A);source=CCDS
```

## How can I let TransVar output sequence context?

The option `--aacontext 5` output +/- 5bp protein sequence context.
```bash
$ transvar ganno -i 'chr17:7577124' --ccds --aacontext 5
```
```text
chr17:7577124	CCDS11118 (protein_coding)	TP53	-
   chr17:g.7577124C>/c.814G>/p.V272	inside_[cds_in_exon_7]
   is_gene_body;aacontext=RNSFE[V]RVCAC;codon_pos=7577122-7577123-7577124;source
   =CCDS
chr17:7577124	CCDS45605 (protein_coding)	TP53	-
   chr17:g.7577124C>/c.814G>/p.V272	inside_[cds_in_exon_7]
   is_gene_body;aacontext=RNSFE[V]RVCAC;codon_pos=7577122-7577123-7577124;source
   =CCDS
chr17:7577124	CCDS45606 (protein_coding)	TP53	-
   chr17:g.7577124C>/c.814G>/p.V272	inside_[cds_in_exon_7]
   is_gene_body;aacontext=RNSFE[V]RVCAC;codon_pos=7577122-7577123-7577124;source
   =CCDS
```
shows the protein sequence context in the aacontext tag.

## How to report results in one line for each query?

Use `--oneline` option. This separates the outputs from each transcript by '|||'.

## I got 'gene_not_recognized', what's wrong?

Most likely you forgot to specify a transcipt definition such as `--ccds` or `--ensembl`. Sometimes there are non-canonical names for genes, this can be fixed through the `--alias` option and specify an alias table. TransVar comes with alias table from UCSC knownGene.

## Does TransVar support alternative format for MNV such as `c.508_509CC>TT`?

Yes, but only in input. For example, `c.508_509CC>TT`
```bash
$ transvar canno --ccds -i 'A1CF:c.508_509CC>TT'
```
```text
A1CF:c.508_509CC>TT	CCDS7241 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   CSQN=Missense;codon_cDNA=508-509-510;source=CCDS
A1CF:c.508_509CC>TT	CCDS7242 (protein_coding)	A1CF	-
   chr10:g.52595929_52595930delinsAA/c.508_509delinsTT/p.P170L	inside_[cds_in_exon_4]
   CSQN=Missense;codon_cDNA=508-509-510;source=CCDS
```

## Does TransVar support relaxed input without 'g.', 'c.' and 'p.'?

Yes, the 'g.', 'c.' and 'p.' are optional in the input. For example, `12:109702119insACC` is equally acceptable as `chr12:g.109702119_109702120insACC`. TransVar also accepts '>' in denoting MNV. E.g., `c.113G>TACTAGC` can be used in place of `c.113delGinsTACTAGC`. This is common in some database such as COSMIC.

## When I annotate a variant for protein identifier, why would I end up getting results in another variant type?

TransVar follows in full the HGVS nomenclature while annotating protein level mutation identifiers. For example, a out-of-phase, in frame insertion, `ACIN1:c.1930_1931insATTCAC` will be annotated with `p.S643_R644insHS` rather than `R644delinsHSR`. Protein level mutation will be generated as if no nucleotide mutation information exists.
