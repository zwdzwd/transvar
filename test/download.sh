## annotation setup ##
transvar config --download_anno --refversion hg38
transvar config --download_anno --refversion hg19
transvar config --download_anno --refversion mm10
transvar config -k reference -v /primary/vari/genomicdata/genomes/mm10/mm10.fa --refversion mm10
transvar config -k reference -v /primary/vari/genomicdata/genomes/hg19/hg19.fa --refversion hg19
transvar config -k refversion -v hg19
transvar config --download_dbsnp
transvar config --download_idmap
