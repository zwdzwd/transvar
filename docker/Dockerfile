FROM python:3.7
MAINTAINER Wanding Zhou <zhouwanding@gmail.com>

ENV TRANSVAR_VERSION==2.4.5.20190121
RUN pip install transvar && \
	transvar config --download_anno --refversion hg38

CMD ['/bin/bash']

