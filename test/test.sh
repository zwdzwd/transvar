################################
## build and upload
################################
# procedures to follow to update

## 1. test
cd test/
python test.py ../docs/source/ outdocs
mv outdocs/*.rst ../docs/source/

## 2. upload github
version="2.5.0.20190314"
# modify ../docs/source/conf.py
# modify ../transvar/version.py
# modify ../docker/Dockerfile

## manually draft a release

## register testpypi
## python setup.py register -r https://testpypi.python.org/pypi
## python setup.py sdist upload -r https://testpypi.python.org/pypi
## pip install -i https://testpypi.python.org/pypi transvar

## 3. upload pypi
## python setup.py register -r https://pypi.python.org/pypi
## python setup.py sdist upload -r https://pypi.python.org/pypi
## pip install -i https://pypi.python.org/pypi transvar

## UPDATE, the current repo is https://upload.pypi.org/legacy/
## /primary/home/wanding.zhou/.local/bin/twine upload dist/*

###### ON HPC environment ####
python setup.py sdist # build the .tar.gz file
# python setup.py bdist # build the .tar.gz file binary (will be platform-specific)
# twine can be found at /primary/home/wanding.zhou/.local/bin/
# or installed from python3 -m pip install --user --upgrade twine
twine upload --repository-url https://upload.pypi.org/legacy/ dist/TransVar-TRANSVARVERSION.tar.gz

# you may want to consider building manually
# docker build -t zhouwanding/transvar:2.4.7 docker

git commit -am "version $version"
git tag -a v$version -m "version $version"
git push
git push --tags

## 4. Update release note on github

# ~/.pypirc
# [distutils]
# index-servers =
#   pypi
#   pypitest

# [pypi]
# repository=https://pypi.python.org/pypi
# username=name
# password=pass

# [pypitest]
# repository=https://testpypi.python.org/pypi
# username=name
# password=pass

############################
## annotation setup for test
############################
transvar config --download_anno --refversion hg38
transvar config --download_anno --refversion hg19
transvar config --download_anno --refversion mm10
transvar config -k reference -v /primary/vari/genomicdata/genomes/mm10/mm10.fa --refversion mm10
transvar config -k reference -v /primary/vari/genomicdata/genomes/hg19/hg19.fa --refversion hg19
transvar config -k refversion -v hg19
transvar config --download_dbsnp
transvar config --download_idmap

####################
## additional test
####################
mkdir -p testout

transvar panno -l cosmic/cosmic_panno_frameshift -g 1 -m 4 --ccds | tee testout/panno_cosmic_frameshift
colordiff testout/
transvar panno -l data/cosmic/cosmic_panno_deletion -g 1 -m 4 --ccds | tee testout/panno_cosmic_deletion
transvar panno -l data/cosmic/cosmic_panno_insertion -g 1 -m 4 --ccds | tee testout/panno_cosmic_insertion
transvar panno -l data/cosmic/cosmic_panno_saav -g 1 -m 4 --ccds | tee testout/panno_cosmic_saav
transvar canno -l data/cosmic/cosmic_canno_snv -g 1 -m 3 --ccds | tee testout/canno_cosmic_snv
transvar canno -l data/cosmic/cosmic_canno_deletion -g 1 -m 3 --ccds | tee testout/canno_cosmic_deletion
transvar canno -l data/cosmic/cosmic_canno_insertion -g 1 -m 3 --ccds | tee testout/canno_cosmic_insertion
colordiff testout/canno_cosmic_insertion golden/canno_cosmic_insertion
transvar panno -l data/tamborero_data/transvar_revprotein_input.txt_ens --ensembl --seqmax -1 | tee testout/panno_tamborero
colordiff testout/panno_tamborero golden/panno_tamborero
transvar ganno --vcf data/tamborero_data/docm_variants_mar_2016.vcf --ensembl | tee testout/ganno_tamborero.vcf
colordiff testout/ganno_tamborero.vcf
transvar ganno -l data/tamborero_data/transvar_dna_input.txt --ensembl | tee testout/ganno_tamborero_output

#####################
## index annotation #
#####################

transvar config --download_raw # check the code for the actual transcript sourcd
transvar index --ccds ~/transvar/Download/hg19.ccds.txt

##################################
## building feature database #####
##################################

### Feature database: Example, hg19 ensembl regulation
cd ~/Dropbox/Public/annotations/hg19/Ensembl75
zcat AnnotatedFeatures.gff.gz | gawk -F"\t" -v OFS="\t" '{match($9,"Cell_type=([^;]*)",a); match($9,"Name=([^;]*)",b); match($9,"Summit=([0-9]*)",c); print $1,$4,$5,b[1]"|"a[1]"|"c[1];}' | sortbed >hg19_AnnotatedFeatures_Ensembl75
transvar index --sorted --bed hg19_AnnotatedFeatures_Ensembl75
rm -f hg19_AnnotatedFeatures_Ensembl75
transvar config -k Histone -v `pwd`/Ensembl75/hg19_AnnotatedFeatures_Ensembl75.featuredb
### Feature database: Example, hg19 Motif
zcat MotifFeatures.gff.gz | awk '{match($9,/Name=([^;]*)/,a); gsub(":","|",a[1]); print $1,$4,$5,"TF_binding_sites|"a[1];}' >hg19_MotifFeatures_Ensembl75
transvar index --bed hg19_MotifFeatures_Ensembl75
rm -f hg19_MotifFeatures_Ensembl75
transvar config -k TFBS -v `pwd`/Ensembl75/hg19_MotifFeatures_Ensembl75.featuredb

### Feature database: Example, hg19 dbsnp
transvar config -k dbsnp -v /Users/wandingzhou/Dropbox/Public/annotations/hg19/hg19_dbsnp.vcf.gz.featuredb

### Feature database: Example, hg38 ensembl regulation
cd ~/Dropbox/Public/annotations/hg38/Ensembl85
zcat AnnotatedFeatures.gff.gz | gawk -F"\t" -v OFS="\t" '{match($9,"epigenome=([^;]*)",a); gsub(" ","_",a[1]); print $1,$4,$5,$3"|"a[1];}' | sortbed >hg38_AnnotatedFeatures_Ensembl85
transvar index --sorted --bed hg38_AnnotatedFeatures_Ensembl85

### Feature database: Example, hg38 Motif
zcat MotifFeatures.gff.gz | awk '{match($9,/motif_feature_type=([^;]*)/,a); match($9, /binding_matrix=([^;]*)/,b); print $1,$4,$5,"TF_binding_sites|"a[1]"|"b[1];}' >hg38_MotifFeatures_Ensembl85
transvar index --bed hg38_MotifFeatures_Ensembl85



