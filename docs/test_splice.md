## test splice site annotation

Consider a splice donor site chr7:5568790_5568791 (a donor site, intron side by definition, reverse strand, chr7:5568792- is the exon),

The 1st exonic nucleotide before donor splice site:
```Shell
$ transvar ganno -i 'chr7:5568792C>G' --ccds
```
output a exonic variation and a missense variation
```text
chr7:5568792C>G	CCDS5341 (protein_coding)	ACTB	-
   chr7:g.5568792C>G/c.363G>C/p.Q121H	inside_[cds_in_exon_2]
   CSQN=Missense;C2=NextToSpliceDonorOfExon2_At_chr7:5568791;codon_pos=5568792-5
   568793-5568794;ref_codon_seq=CAG;source=CCDS
```

The 1st nucleotide in the canonical donor splice site (intron side, this is commonly regarded as the splice site location):
```Shell
$ transvar ganno -i 'chr7:5568791C>G' --ccds
```
output a splice variation
```text
chr7:5568791C>G	CCDS5341 (protein_coding)	ACTB	-
   chr7:g.5568791C>G/c.363+1G>C/.	inside_[intron_between_exon_2_and_3]
   CSQN=SpliceDonorSNV;C2=SpliceDonorOfExon2_At_chr7:5568791;source=CCDS
```

The 2nd nucleotide in the canonical donor splice site (2nd on the intron side, still considered part of the splice site):
```Shell
$ transvar ganno -i 'chr7:5568790A>G' --ccds
```
output a splice variation
```text
chr7:5568790A>G	CCDS5341 (protein_coding)	ACTB	-
   chr7:g.5568790A>G/c.363+2T>C/.	inside_[intron_between_exon_2_and_3]
   CSQN=SpliceDonorSNV;C2=SpliceDonorOfExon2_At_chr7:5568791;source=CCDS
```

The 1st nucleotide downstream next to the canonical donor splice site (3rd nucleotide in the intron side, not part of the splice site):
```Shell
$ transvar ganno -i 'chr7:5568789C>G' --ccds
```
output a pure intronic variation
```text
chr7:5568789C>G	CCDS5341 (protein_coding)	ACTB	-
   chr7:g.5568789C>G/c.363+3G>C/.	inside_[intron_between_exon_2_and_3]
   CSQN=IntronicSNV;source=CCDS
```

