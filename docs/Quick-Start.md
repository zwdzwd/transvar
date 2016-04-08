
Here we show how one can use TransVar on human hg19 (GRCh37). 
```bash
# set up databases
transvar config --download_anno --refversion hg19
# in case you don't have a reference
transvar config --download_ref --refversion hg19
# in case you do have a reference to link
transvar config -k reference -v [path_to_hg19.fa] --refversion hg19
```

Test an input:
```bash
$ transvar panno -i 'PIK3CA:p.E545K' --ucsc --ccds
```
outputs show two hits from the two databases, i.e., UCSC and CCDS.
```text
PIK3CA:p.E545K	NM_006218 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	inside_[cds_in_exon_10]
   CSQN=Missense;reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_vari
   ants=chr3:g.178936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G
   >A);source=UCSCRefGene
PIK3CA:p.E545K	CCDS43171 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	inside_[cds_in_exon_9]
   CSQN=Missense;reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_vari
   ants=chr3:g.178936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G
   >A);source=CCDS
```
One could provide input based on transcript ID, e.g `NM_006218.1:p.E545K` and TransVar would automatically restrict to the provided transcript.
```
$ transvar panno -i 'NM_006218.2:p.E545K' --ucsc --ccds
```
outputs
```text
NM_006218.2:p.E545K	NM_006218 (protein_coding)	PIK3CA	+
   chr3:g.178936091G>A/c.1633G>A/p.E545K	inside_[cds_in_exon_10]
   CSQN=Missense;reference_codon=GAG;candidate_codons=AAG,AAA;candidate_mnv_vari
   ants=chr3:g.178936091_178936093delGAGinsAAA;dbsnp=rs104886003(chr3:178936091G
   >A);source=UCSCRefGene
```

