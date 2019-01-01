
check dbsnp
```Shell
$ transvar ganno -i 'chr1:54721del' --ccds
```
```text
chr1:54721del	.	.	.
   chr1:g.54723delT/./.	inside_[intergenic_between_5'-telomere(54,721_bp)_and_OR4F5(14,370_bp_upstream)]
   CSQN=IntergenicDeletion;dbsnp=rs750488165(chr1:54720CTTTCTTTCTTTCTTTCT>C),rs3
   73430935(chr1:54720CT>C),rs765718081(chr1:54720CTTTCT>C);left_align_gDNA=g.54
   721delT;unaligned_gDNA=g.54721delT
```

check dbsnp
```Shell
$ transvar ganno -i 'chr1:1066952_1066953AT>GC' --ccds
```
```text
chr1:1066952_1066953AT>GC	.	.	.
   chr1:g.1066952_1066953delinsGC/./.	inside_[intergenic_between_C1ORF159(40,029_bp_upstream)_and_TTLL10(47,643_bp_upstream)]
   CSQN=IntergenicBlockSubstitution;dbsnp=rs34955020(chr1:1066951CAT>CGC)
```

deletion + snv
```Shell
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCTACCCAGGAG' --haplotype
```
```text
chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCTACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645101delG/c.138delC/p.Y46*fs*1	inside_[cds_in_exon_2]
   CSQN=Frameshift;left_align_gDNA=g.645101delG;unaligned_gDNA=g.645101delG;left
   _align_cDNA=c.138delC;unalign_cDNA=c.138delC;source=CCDS
chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCTACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	inside_[intron_between_exon_1_and_2]
   CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon1_At_chr20:645106;source=CCDS
```

insertion + snv
```Shell
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCGATACCCAGGAG' --haplotype
```
```text
chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350	inside_[cds_in_exon_2]
   CSQN=Frameshift;left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645
   102insA;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;source=CCDS
chr20:g.645097_645111GGGCGTACCCTGGAG>GGGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	inside_[intron_between_exon_1_and_2]
   CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon1_At_chr20:645106;source=CCDS
```

snv + snv
```Shell
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGTACCCAGGAG' --haplotype
```
```text
chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGTACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645098G>T/c.141C>A/p.A47A	inside_[cds_in_exon_2]
   CSQN=Synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC;source=CCDS
chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGTACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	inside_[intron_between_exon_1_and_2]
   CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon1_At_chr20:645106;source=CCDS
```

snv + insertion + snv
```Shell
$ transvar ganno --ccds -i 'chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGATACCCAGGAG' --haplotype
```
```text
chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645098G>T/c.141C>A/p.A47A	inside_[cds_in_exon_2]
   CSQN=Synonymous;codon_pos=645098-645099-645100;ref_codon_seq=GCC;source=CCDS
chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645101_645102insA/c.137_138insT/p.A47Rfs*350	inside_[cds_in_exon_2]
   CSQN=Frameshift;left_align_gDNA=g.645101_645102insA;unalign_gDNA=g.645101_645
   102insA;left_align_cDNA=c.137_138insT;unalign_cDNA=c.137_138insT;source=CCDS
chr20:g.645097_645111GGGCGTACCCTGGAG>GTGCGATACCCAGGAG	CCDS13006 (protein_coding)	SCRT2	-
   chr20:g.645107T>A/c.134-2A>T/.	inside_[intron_between_exon_1_and_2]
   CSQN=SpliceAcceptorSNV;C2=SpliceAcceptorOfExon1_At_chr20:645106;source=CCDS
```

### test reduction format relaxation

Frameshift without first alternative allele or termination length
```Shell
$ transvar panno -i 'APC:p.I1557fs' --ccds
```
```text
APC:p.I1557fs	CCDS4107 (protein_coding)	APC	+
   chr5:g.(112175959_112175960)/c.(4669_4668)/p.I1557fs	inside_[cds_in_exon_15]
   CSQN=Frameshift;imprecise;source=CCDS
```

Frameshift with first alternative allele (a *) but no termination length
```Shell
$ transvar panno -i 'APC:p.I1557*fs' --ccds
```
```text
APC:p.I1557*fs	CCDS4107 (protein_coding)	APC	+
   chr5:g.112175960_112175962delATTinsTAA/c.4669_4671delATTinsTAA/p.I1557*	inside_[cds_in_exon_15]
   CSQN=Nonsense;reference_codon=ATT;candidate_codons=TAA,TAG,TGA;candidate_mnv_
   variants=chr5:g.112175960_112175962delATTinsTAG,chr5:g.112175960_112175962del
   ATTinsTGA;source=CCDS
```

Frameshift with termination length but no first alternative allele 
```Shell
$ transvar panno -i 'APC:p.I1557fs*3' --ccds
```
```text
APC:p.I1557fs*3	CCDS4107 (protein_coding)	APC	+
   chr5:g.112175961_112175962delTT/c.4670_4671delTT/p.I1557fs*3	inside_[cds_in_exon_15]
   CSQN=Frameshift;left_align_cDNA=c.4670_4671delTT;left_align_gDNA=g.112175961_
   112175962delTT;source=CCDS
```

Frameshift with first alternative allele but no termination length
```Shell
$ transvar panno -i 'APC:p.I1557Efs' --ccds
```
```text
APC:p.I1557Efs	CCDS4107 (protein_coding)	APC	+
   chr5:g.(112175959_112175960)/c.(4669_4668)/p.I1557Efs	inside_[cds_in_exon_15]
   CSQN=Frameshift;imprecise;source=CCDS
```

```Shell
transvar ganno --ccds -i 'chr3:g.178936091G>A' --gseq
