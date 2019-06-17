FROM python:3.7
MAINTAINER Wanding Zhou <zhouwanding@gmail.com>

ENV TRANSVAR_VERSION==2.5.7.20190617
ENV TRANSVAR_CFG=/anno/transvar.cfg
ENV TRANSVAR_DOWNLOAD_DIR=/anno/
RUN cd / && \
    wget "https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2" && \
    tar -jxvf samtools-1.3.1.tar.bz2 && \
    cd samtools-1.3.1 && \
    make && \
    cp samtools /usr/bin && \
    cd / && \
    wget "https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2" && \
    tar -jxvf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make && \
    cp bgzip tabix /usr/bin/ && \
    cd / && \
	pip install transvar==${TRANSVAR_VERSION} && \
	transvar config --download_anno --refversion hg38 --skip_reference

VOLUME [ "/data" ]

CMD ["transvar"]

## this how to run this
## assuming the existence of ~/references/hg38/hg38.fa and ~/references/hg38/hg38.fa.fai
## docker run -v ~/references/hg38:/data -ti zhouwanding/transvar transvar panno -i PIK3CA:p.E545K --ensembl --reference /data/hg38.fa
