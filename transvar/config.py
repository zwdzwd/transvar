"""
The MIT License

Copyright (c) 2015
The University of Texas MD Anderson Cancer Center
Wanding Zhou, Tenghui Chen, Ken Chen (kchen3@mdanderson.org)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import configparser
import subprocess
import os, sys
from .err import *

from future.standard_library import install_aliases
install_aliases()
from urllib.request import urlopen
# import urllib.request, urllib.error, urllib.parse

samtools_path='%s/samtools' % os.path.abspath(os.path.dirname(__file__))

def samtools_faidx(fn):
    
    err_print("Faidx indexing")
    subprocess.check_call([samtools_path, 'faidx', fn])

def gunzip(fn):

    if not fn.endswith('.gz'):
        err_die('Target file %s not ends with .gz' % fn)

    import gzip
    f_out = open(fn[:-3], 'w')
    f_in = gzip.open(fn)
    f_out.writelines(f_in)
    f_in.close()
    f_out.close()
    os.remove(fn)

    
cfg_fns = [
    os.path.expanduser(os.getenv('TRANSVAR_CFG', 
                                 os.path.join(os.path.dirname(__file__), 'transvar.cfg'))),
    os.path.expanduser('~/.transvar.cfg')]

downloaddirs = [
    os.path.expanduser(os.getenv('TRANSVAR_DOWNLOAD_DIR',
                                 os.path.join(os.path.dirname(__file__), 'transvar.download'))),
    os.path.expanduser('~/.transvar.download')]

# dwroot = 'https://dl.dropboxusercontent.com/u/6647241/annotations/'
dwroot = 'http://transvar.info/transvar_user/annotations/'

fns = {}

# download link updates:
# 03/06/2016: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz
fns[('hg19', 'raw')] = [
    ('refseq', 'hg19.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz'),
    ('ccds', 'hg19.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt'),
    ('ensembl', 'hg19.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'),
    ('gencode', 'hg19.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'),
    ('ucsc', 'hg19.ucsc.txt.gz', '%s/hg19.ucsc.refgene.txt.gz' % dwroot),
    # ('custom', 'hg19.custom.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.map?dl=1'),
    ('aceview', 'hg19.aceview.gff.gz', 'ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Aug10.human.genes/AceView.ncbi_37.genes_gff.gff.gz'),
    ('known_gene', 'hg19.knowngene.gz', '%s/UCSC_knownGene_hg19.gz?dl=1' % dwroot),
    (None, 'hg19.knowngene_alias.gz', '%s/UCSC_kgAlias.gz?dl=1' % dwroot),
]

fns[('hg19', 'dbsnp')] = [
    ('dbsnp', 'hg19_dbsnp.vcf.gz', '%s/hg19_dbsnp.vcf.gz' % dwroot),
    (None, 'hg19_dbsnp.vcf.gz.tbi', '%s/hg19_dbsnp.vcf.gz.tbi' % dwroot),
    # 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/00-All.vcf.gz'),
    # (None, 'hg19_dbsnp.vcf.gz.tbi', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/00-All.vcf.gz.tbi'),
]

fns[('hg18', 'raw')] = [
    ('refseq', 'hg18.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/GFF/ref_NCBI36_top_level.gff3.gz'),
    ('ccds', 'hg18.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs36.3/CCDS.20090327.txt'),
    # the AceView hg18 version is deprecated.
    # ('aceview', 'hg18.aceview.gff.gz', 'ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_36_Apr07.human.genes/AceView.ncbi_36.genes_gff.tar.gz'),
    ('gencode', 'hg18.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_3c/gencode.v3c.annotation.NCBI36.gtf.gz'),
    ('ucsc', 'hg18.ucsc.txt.gz', '%s/hg18.ucsc.refgene.txt.gz?dl=1' % dwroot),
    ('ensembl', 'hg18.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz'),
]

for rv in ['hg18', 'hg19', 'hg38', 'mm9', 'mm10']:
    fns[(rv,'reference')] = [
        ('reference', '%s.fa' % rv, '%s/%s.fa' % (dwroot, rv)),
        (None, '%s.fa.fai' % rv, '%s/%s.fa.fai' % (dwroot, rv)),
    ]

# download link update:
# 01/05/2015: refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p2_top_level.gff3.gz
# 01/05/2015: ensembl: ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
# 06/27/2016: refseq: ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz
# 06/27/2016: ensembl: ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
# 06/27/2016: gencode: ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
fns[('hg38', 'raw')] = [
    ('refseq', 'hg38.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p7_top_level.gff3.gz'),
    ('ccds', 'hg38.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt'),
    ('ensembl', 'hg38.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz'),
    ('gencode', 'hg38.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz'),
    ('ucsc', 'hg38.ucsc.txt.gz', '%s/hg38.ucsc.refgene.txt.gz?dl=1' % dwroot),
]

fns[('mm9', 'raw')] = [
    ('ensembl', 'mm9.ensembl.gtf.gz',
     'ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz'),
    ('ccds', 'mm9.ccds.txt',
     'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Mm37.1/CCDS.current.txt'),
    ('gencode', 'mm9.gencode.gtf.gz',
     'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz'),
]

fns[('mm10', 'raw')] = [
    ('refseq', 'mm10.refseq.gff.gz',
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/GFF/ref_GRCm38.p3_top_level.gff3.gz'),
    ('ccds', 'mm10.ccds.txt',
     'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Mm38.1/CCDS.current.txt'),
    ('ensembl', 'mm10.ensembl.gtf.gz',
     'ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz'),
    ('gencode', 'mm10.gencode.gtf.gz',
     'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz'), # GRCm38.p3 genome
]

# build anno topics
fns2 = []
for (refv, topic), vs in fns.items():
    if topic == 'raw':          # convert all the raw files
        vs2 = []
        for k, fn, url in vs:
            if k is None:
                continue
            afn = fn+'.transvardb'
            vs2.append((k, afn, dwroot+afn))
            afn = fn+'.transvardb.gene_idx'
            vs2.append((None, afn, dwroot+afn))
            afn = fn+'.transvardb.trxn_idx'
            vs2.append((None, afn, dwroot+afn))
            afn = fn+'.transvardb.loc_idx'
            vs2.append((None, afn, dwroot+afn))
            afn = fn+'.transvardb.loc_idx.tbi'
            vs2.append((None, afn, dwroot+afn))
            # if (k.find('knowngene')>=0 or k.find('refseq')>=0 or k.find('gencode')>=0 or k.find('ensembl')>=0):
            afn = fn+'.transvardb.alias_idx'
            vs2.append((None, afn, dwroot+afn))

        fns2.append(((refv, 'anno'), vs2))

for k, v in fns2:
    fns[k] = v

# def download_hg19_reference(config):
#     fns = [('reference', 'hg19.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)

# def download_hg38_reference(config):
#     fns = [('reference', 'hg38.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)

def download_url(url, file_name):

    import ssl
    # file_name = url.split('/')[-1]
    # try:
    if hasattr(ssl, '_create_unverified_context'):
        u = urlopen(url, context=ssl._create_unverified_context())
    else:
        u = urlopen(url)

    # except urllib2.URLError:
    # return
    f = open(file_name, 'wb')
    meta = u.info()
    raw_file_size = int(meta.getheaders("Content-Length")[0])
    file_size = raw_file_size / (1024.0 * 1024.0)
    # err_print("downloading %s (%1.1f MB)" % (file_name, file_size))

    file_size_dl = 0
    block_sz = 8192*2
    sys.stdout.write('[downloading] %s ...' % (file_name, ))
    sys.stdout.flush()
        
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        # status = r"downloaded %s (%1.1f MB) %10d [%3.2f%%]\033\[K" % (file_name, file_size, file_size_dl, file_size_dl * 100. / raw_file_size)
        # status = status + chr(8)*(len(status)+1)
        progress = float(file_size_dl)/raw_file_size

    print('Done (%1.1f MB).' % (file_size_dl/1000000., ))
    f.close()

def download_requests(url, file_name):

    """ sometimes, urllib2 doesn't work, try requests """
    import requests
    r = requests.get(url, stream=True)
    if r.status_code != 404:
        sys.stdout.write('[bakdownloading] %s ...' % (file_name, ))
        sys.stdout.flush()
        with open(file_name,'wb') as fd:
            n = 0
            
            for chunk in r.iter_content(10000000):
                n += len(chunk)
                fd.write(chunk)

        print('Done (%1.1f MB).' % (n/1000000., ))

def config_set(config, section, option, value):

    if section != 'DEFAULT' and not config.has_section(section):
        config.add_section(section)
    config.set(section, option, value)

def _download_(config, section, fns):

    for pdir in downloaddirs:
        # pdir = os.path.join(os.path.dirname(__file__), 'download')

        if not os.path.exists(pdir):
            try:
                os.makedirs(pdir)
            except:
                continue

        for k, fn, url in fns:
            fnn = os.path.join(pdir, fn)
            
            success = True
            try:
                download_url(url, fnn)
                if k:
                    config_set(config, section, k, fnn)
            except:
                if not fn.endswith('alias_idx'): # sometimes some alias_idx will be missing
                    success = False

            if success:
                continue

            try:                # not quite necessary in most cases, but in some situations, urllib won't work
                download_requests(url, fnn)
                if k:
                    config_set(config, section, k, fnn)
            except:
                err_warn('file not available: %s or target directory not found' % url)

        break
    return pdir

def download_idmap(config):
    # 'https://dl.dropboxusercontent.com/u/6647241/annotations/HUMAN_9606_idmapping.dat.gz?dl=1'
    fns = [('uniprot', 'uniprot.idmapping.txt.gz.idx',
            '%s/uniprot.idmapping.txt.gz.idx' % dwroot)]
    _download_(config, 'idmap', fns)

def getrv(args, config):

    if args.refversion != 'DEFAULT':
        rv = args.refversion
    elif 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
    else:
        rv = 'hg19'

    return rv

def download_topic(args, config, topic):

    rv = getrv(args, config)
    if (rv, topic) in fns:
        config.set('DEFAULT', 'refversion', rv)
        _download_(config, rv, fns[(rv, topic)])
    else:
        err_die('no pre-built %s for %s, please build manually' % (topic, rv))


def download_anno_topic_ensembl(args, config):

    from ftplib import FTP
    rv = getrv(args, config)
    args.ensembl_release = 80
    eshost = 'ftp.ensembl.org'
    ftp = FTP(eshost)
    ftp.login()

    esroot = 'pub/release-%d/' % args.ensembl_release
    if args.refversion == 'DEFAULT':
        species = [os.path.basename(o) for o in ftp.nlst("%s/gtf/" % esroot)]
        for i, sp in enumerate(species):
            err_print('[%d] %s' % (i,sp))
        choice = input("Please choose your target taxon [%d]: " % species.index("homo_sapiens"))
        if not choice:
            choice = species.index('homo_sapiens')
        choice = int(choice)
        err_print("Preparing genomes and annotations for [%d] %s." % (choice, species[choice]))
        if int(choice) <= 0 or int(choice) >= len(species):
            err_die("Invalid choice.")
        rv = species[choice]
    
    esfasta = '%s/fasta/%s/dna/' % (esroot, rv.lower())
    genomes = [fn for fn in ftp.nlst(esfasta) if fn.endswith('dna.toplevel.fa.gz')]
    assert(len(genomes) == 1)
    genome = genomes[0]
    genoname = os.path.basename(genome)
    genodir = _download_(config, rv, [(None, genoname, 'ftp://'+eshost+'/'+genome)])

    esgtf = '%s/gtf/%s' % (esroot, rv.lower())
    gtfs = [fn for fn in ftp.nlst(esgtf) if fn.endswith('gtf.gz')]
    assert(len(gtfs)==1)
    gtf = gtfs[0]
    gtfname = os.path.basename(gtf)
    gtfdir = _download_(config, rv, [(None, gtfname, 'ftp://'+eshost+'/'+gtf)])
    
    err_print("Unzipping genome")
    gunzip(genodir+'/'+genoname)

    samtools_faidx(genodir+'/'+genoname[:-3])
    config_set(config, rv, 'reference', genodir+'/'+genoname[:-3])

    err_print("Indexing GTF")
    from . import localdb
    db = localdb.EnsemblDB()
    db.index([gtfdir+'/'+gtfname])
    config_set(config, rv, 'ensembl', gtfdir+'/'+gtfname+'.transvardb')

    config.set('DEFAULT', 'refversion', rv)

def read_config():
    config = configparser.RawConfigParser()
    config.read(cfg_fns)
    return config

def main(args):

    config = configparser.RawConfigParser()
    config.read(cfg_fns)

    if args.k and args.v:
        if args.k == 'refversion':
            sec = 'DEFAULT'
        else:
            sec = getrv(args, config)
        config_set(config, sec, args.k, args.v)
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)

    if args.download_ref:
        download_topic(args, config, 'reference')
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)

    if args.download_anno:
        download_topic(args, config, 'anno')
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)

    if args.download_ensembl:
        download_anno_topic_ensembl(args, config)
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)

    if args.download_raw:
        download_topic(args, config, 'raw')
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)
            
    if args.download_dbsnp:
        download_topic(args, config, 'dbsnp')
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)
        
    if args.download_idmap:
        download_idmap(config)
        if args.refversion != 'DEFAULT':
            config.set('DEFAULT', 'refversion', args.refversion)

    for cfg_fn in cfg_fns:
        try:
            config.write(open(cfg_fn,'w'))
            break
        except IOError as e:
            pass

def main_current(args):

    config = configparser.RawConfigParser()
    config.read(cfg_fns)
    if args.refversion != 'DEFAULT':
        rv = args.refversion
    elif 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
    else:
        err_die("no default reference version set.")

    print("reference version: %s" % rv)
    if 'reference' in config.options(rv):
        print('reference: %s' % config.get(rv, 'reference'))
    print("Available databases: ")
    for op in config.options(rv):
        if op not in ['refversion', 'reference']:
            print('%s: %s' % (op, config.get(rv, op)))

    return


def add_parser_config(subparsers):

    parser = subparsers.add_parser('config', help="show configurations")
    parser.add_argument('-k', default=None, help='key')
    parser.add_argument('-v', default=None, help='set value')
    parser.add_argument('--refversion', default='DEFAULT',
                        help='reference version, options: hg18, hg19, hg38, mm9, mm10, see transvar config --download_ensembl for others.')
    parser.add_argument('--download_anno', action='store_true', help='download annotations')
    parser.add_argument('--download_ensembl', action='store_true', help='download ensembl raw annotations')
    parser.add_argument('--download_ref', action='store_true', help='download reference')
    parser.add_argument('--download_dbsnp', action='store_true', help='download dbsnp')
    parser.add_argument('--download_idmap', action='store_true', help='download id map')
    parser.add_argument('--download_raw', action='store_true', help='download annotation raw file')
    parser.set_defaults(func=main)

def add_parser_current(subparsers):

    parser = subparsers.add_parser('current', help="view current config")
    parser.add_argument('--refversion', default='DEFAULT', help='reference version')
    parser.set_defaults(func=main_current)
