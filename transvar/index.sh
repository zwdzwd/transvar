
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

for rv in hg19 hg18 hg38 mm10 mm9; do
  transvar config --download_anno --refversion $rv;
done
