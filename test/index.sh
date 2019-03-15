
# Note for indexing
for rv in hg19 hg38 mm10; do
  transvar config --download_raw --refversion $rv;
done

cd transvar.download
for f in *;do
  if [[ $f != *transvardb* ]]; then
    echo $f;
    if [[ $f == *ccds* ]]; then
      transvar index --ccds $f;
    fi
    if [[ $f == *ensembl* ]]; then
      transvar index --ensembl $f;
    fi
    if [[ $f == *gencode* ]]; then
      transvar index --gencode $f;
    fi
    if [[ $f == *refseq* ]]; then
      transvar index --refseq $f;
    fi
    if [[ $f == *aceview* ]]; then
      transvar index --aceview $f;
    fi
    if [[ $f == *ucsc* ]]; then
      transvar index --ucsc $f;
    fi
    if [[ $f == *knowngene.gz* ]]; then
      transvar index --kg $f --alias ${f%%.gz}_alias.gz;
    fi
    if [[ $f == *uniprot* ]]; then
      transvar index --uniprot $f;
    fi
  fi
done

for rv in hg19 hg38 mm10; do # mm9 and hg18 are not supported any more
  transvar config --download_ref --refversion $rv;
  transvar config --download_anno --refversion $rv;
done


# hg19.aceview.gff.gz
# [parse_raw] loaded 259440 transcripts from AceView GFF file.
# hg19.ccds.txt
# [parse_raw] loaded 26283 transcripts from CCDS table.
# hg19.ensembl.gtf.gz
# [parse_raw] loaded 215170 transcripts from Ensembl GTF file.
# hg19.gencode.gtf.gz
# [parse_raw] loaded 196520 transcripts from GENCODE GTF file.
# hg19.knowngene.gz
# [parse_raw] loaded 82960 transcripts from UCSC knownGene table.
# hg19.knowngene_alias.gz
# hg19.refseq.gff.gz
# [parse_raw] loaded 68278 transcripts from RefSeq GFF3 file.
# hg19.ucsc.txt.gz
# [parse_raw] loaded 40634 transcripts from UCSC refgene.
# hg38.ccds.txt
# [parse_raw] loaded 33381 transcripts from CCDS table.
# hg38.ensembl.gtf.gz
# [parse_raw] loaded 206601 transcripts from Ensembl GTF file.
# hg38.gencode.gtf.gz
# [parse_raw] loaded 206694 transcripts from GENCODE GTF file.
# hg38.refseq.gff.gz
# [parse_raw] loaded 113775 transcripts from RefSeq GFF3 file.
# hg38.ucsc.txt.gz
# [parse_raw] loaded 43231 transcripts from UCSC refgene.
# mm10.ccds.txt
# [parse_raw] loaded 22934 transcripts from CCDS table.
# mm10.ensembl.gtf.gz
# [parse_raw] loaded 138930 transcripts from Ensembl GTF file.
# mm10.gencode.gtf.gz
# [parse_raw] loaded 138835 transcripts from GENCODE GTF file.
# mm10.refseq.gff.gz
# [parse_raw] loaded 100205 transcripts from RefSeq GFF3 file.
