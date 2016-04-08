TransVar uses optional additional resources for annotation.

## dbSNP

For example, one could annotate SNP with dbSNP id by downloading the dbSNP files.
This can be done by
```bash
transvar config --download_dbsnp
```
TransVar automatically download dbSNP file which correspoding to the current default reference version (as set in `transvar.cfg`). This also sets the entry in `transvar.cfg`.
With dbSNP file downloaded, TransVar automatically looks for dbSNP id when performing annotation.
```bash
$ transvar panno -i 'A1CF:p.A309A' --ccds
```
```text
A1CF:p.A309A	CCDS7243 (protein_coding)	A1CF	-
   chr10:g.52576004T>G/c.927A>C/p.A309A	inside_[cds_in_exon_7]
   CSQN=Synonymous;reference_codon=GCA;candidate_codons=GCC,GCG,GCT;candidate_sn
   v_variants=chr10:g.52576004T>C,chr10:g.52576004T>A;dbsnp=rs201831949(chr10:52
   576004T>G);source=CCDS
```
Note that in order to use dbSNP, one must download the dbSNP database through `transvar config --download_dbsnp`, or by configure the `dbsnp` slot in the configure file via `transvar config -k dbsnp -v [path to dbSNP VCF]`. Manually set path for dbSNP file must have the file tabix indexed.
