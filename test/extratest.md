
check dbsnp
$ transvar ganno -i 'chr1:54721del' --ccds
chr1:54721del	.	.	.
   chr1:g.54723delT/./.	intergenic_between_5'-telomere(54,721_bp)_and_OR4F5(14,370_bp_upstream)
   left_align_gDNA=g.54721delT;unaligned_gDNA=g.54721delT;dbsnp=rs373430935(chr1
   :54720CT>C)

check dbsnp
$ transvar ganno -i 'chr1:1066952_1066953AT>GC' --ccds
chr1:1066952_1066953AT>GC	.	.	.
   chr1:g.1066952_1066953AT>GC/./.	inside_[intergenic_between_C1ORF159(40,029_bp_upstream)_and_TTLL10(47,643_bp_upstream)]
   dbsnp=rs34955020(chr1:1066951CAT>CGC),rs386627447(chr1:1066951CAT>CGC)

deletion + snv
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCTACCCAGGAG' --haplotype
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645101delG/c.138delC/p.Y46*fs*1
   cds_in_exon_2	left_align_gDNA=g.645101delG;unaligned_gDNA=g.645101delG;left_align_cDNA=c.138delC;unalign_cDNA=c.138delC
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645107T>A/c.134-2A>T/.
   intron_between_exon_1_and_2	acceptor_splice_site_of_exon_1_at_chr20:645106

insertion + snv
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCGATACCCAGGAG' --haplotype
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350
   cds_in_exon_2	left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645102insA;insertion_gDNA=A;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;insertion_cDNA=T
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645107T>A/c.134-2A>T/.
   intron_between_exon_1_and_2	acceptor_splice_site_of_exon_1_at_chr20:645106

snv + snv
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGTACCCAGGAG' --haplotype
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645098G>T/c.141C>A/p.A47A
   cds_in_exon_2	synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645107T>A/c.134-2A>T/.
   intron_between_exon_1_and_2	acceptor_splice_site_of_exon_1_at_chr20:645106

snv + insertion + snv
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGATACCCAGGAG' --haplotype
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645098G>T/c.141C>A/p.A47A
   cds_in_exon_2	synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350
   cds_in_exon_2	left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645102insA;insertion_gDNA=A;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;insertion_cDNA=T
CCDS13006.1 (protein_coding)	SCRT2	-	chr20:g.645107T>A/c.134-2A>T/.
   intron_between_exon_1_and_2	acceptor_splice_site_of_exon_1_at_chr20:645106

