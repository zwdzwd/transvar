""" configure TransVar """
import ConfigParser
import os, sys

cfg_fns = [os.path.join(os.path.dirname(__file__), 'transvar.cfg'),
           os.path.expanduser('~/.transvar.cfg')]

downloaddirs = [os.path.join(os.path.dirname(__file__), 'transvar.download'),
                os.path.expanduser('~/.transvar.download')]

fns = {}

fns[('hg19', 'anno')] = [
    ('refseq', 'hg19.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/ANNOTATION_RELEASE.105/GFF/ref_GRCh37.p13_top_level.gff3.gz'),
    ('ccds', 'hg19.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs37.3/CCDS.current.txt'),
    ('ensembl', 'hg19.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz'),
    ('gencode', 'hg19.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz'),
    ('ucsc', 'hg19.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.ucsc.refgene.txt.gz?dl=1'),
    ('custom', 'hg19.custom.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.map?dl=1'),
    ('aceview', 'hg19.aceview.gff.gz', 'ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_37_Aug10.human.genes/AceView.ncbi_37.genes_gff.gff.gz'),
    ('known_gene', 'hg19.knowngene.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_knownGene_hg19.gz?dl=1'),
    ('known_gene_alias', 'hg19.knowgene_alias.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_kgAlias.gz?dl=1'),
]

fns[('hg19', 'dbsnp')] = [
    ('dbsnp', 'hg19_dbsnp.vcf.gz', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz'),
    ('dbsnp_index', 'hg19_dbsnp.vcf.gz.tbi', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz.tbi'),
]

fns[('hg18', 'anno')] = [
    ('refseq', 'hg18.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/GFF/ref_NCBI36_top_level.gff3.gz'),
    ('ccds', 'hg18.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs36.3/CCDS.20090327.txt'),
    ('aceview', 'hg18.aceview.gff.gz', 'ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_36_Apr07.human.genes/AceView.ncbi_36.genes_gff.tar.gz'),
    ('gencode', 'hg18.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_3c/gencode.v3c.annotation.NCBI36.gtf.gz'),
    ('ucsc', 'hg18.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg18.ucsc.refgene.txt.gz?dl=1'),
    ('ensembl', 'hg18.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz'),
]

fns[('hg38', 'anno')] = [
    ('refseq', 'hg38.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38_top_level.gff3.gz'),
    # ('ccds', 'hg38.ccds.txt', ''),
    ('ensembl', 'hg38.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz'),
    ('gencode', 'hg38.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz'),
    ('ucsc', 'hg38.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg38.ucsc.refgene.txt.gz?dl=1'),
]

fns[('mm9', 'anno')] = [
    ('refseq', 'mm9.refseq.gff.gz',
     'ftp://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz'),
    ('ccds', 'mm9.ccds.txt',
     'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Mm37.1/CCDS.current.txt'),
    ('gencode', 'mm9.gencode.gtf.gz',
     'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz'),
]

fns[('mm10', 'anno')] = [
    ('refseq', 'mm10.refseq.gff.gz',
     'ftp://ftp.ncbi.nlm.nih.gov/genomes/M_musculus/GFF/ref_GRCm38.p3_top_level.gff3.gz'),
    ('ccds', 'mm10.ccds.txt',
     'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Mm38.1/CCDS.current.txt'),
    ('ensembl', 'mm10.ensembl.gtf.gz',
     'ftp://ftp.ensembl.org/pub/release-79/gtf/mus_musculus/Mus_musculus.GRCm38.79.gtf.gz'),
    ('gencode', 'mm10.gencode.gtf.gz',
     'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M4/gencode.vM4.annotation.gtf.gz'), # GRCm38.p3 genome
]

# def download_hg19_reference(config):
#     fns = [('reference', 'hg19.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)
    
# def download_hg38_reference(config):
#     fns = [('reference', 'hg38.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)

def download_url(url, file_name):

    import urllib2
    # file_name = url.split('/')[-1]
    u = urllib2.urlopen(url)
    f = open(file_name, 'wb')
    meta = u.info()
    raw_file_size = int(meta.getheaders("Content-Length")[0])
    file_size = raw_file_size / (1024.0 * 1024.0)
    print "Downloading: %s (%1.1f MB)" % (file_name, file_size)

    file_size_dl = 0
    block_sz = 8192*2
    while True:
        buffer = u.read(block_sz)
        if not buffer:
            break

        file_size_dl += len(buffer)
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / raw_file_size)
        status = status + chr(8)*(len(status)+1)
        print status,

    f.close()

def config_set(config, section, option, value):

    if section != 'DEFAULT' and not config.has_section(section):
        config.add_section(section)
    config.set(section, option, value)

def _download_(config, section, fns):

    for pdir in downloaddirs:
        # pdir = os.path.join(os.path.dirname(__file__), 'download')

        try:
            if not os.path.exists(pdir): os.makedirs(pdir)
            for k, fn, url in fns:
                fnn = os.path.join(pdir, fn)
                download_url(url, fnn)
                if k:
                    config_set(config, section, k, fnn)
        except:
            continue

        break

def download_idmap(config):
    fns = [('uniprot', 'uniprot.idmapping.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/HUMAN_9606_idmapping.dat.gz?dl=1')]
    _download_(config, 'idmap', fns)

def getrv(config):

    if args.refversion:
        rv = args.refversion
    elif args.s != 'DEFAULT':
        rv = args.s
    elif 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
    else:
        rv = 'hg19'

    return rv

def download_topic(config, topic):

    rv = getrv(config)
    if (rv, topic) in fns:
        config.set('DEFAULT', 'refversion', fns[(rv, topic)])
        _download_(config, rv, fns[(rv, topic)])
    else:
        err_die('no pre-built %s for %s, please build manually' % (topic, rv))

def read_config():
    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    return config

def main(args):

    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    if args.k and args.v:
        config_set(config, args.s, args.k, args.v)

    if args.download_anno:
        download_topic(config, 'anno')

    if args.download_dbsnp:
        download_topic(config, 'dbsnp')    
        
    if args.download_idmap:
        download_idmap(config)

    for cfg_fn in cfg_fns:
        try:
            config.write(open(cfg_fn,'w'))
            break
        except IOError as e:
            pass

def main_current(args):

    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    if 'refversion' in config.defaults():
        rv = config.get('DEFAULT', 'refversion')
        print "Current reference version: %s" % rv
        if 'reference' in config.options(rv):
            print 'reference: %s' % config.get(rv, 'reference')
        print "Available databases: "
        for op in config.options(rv):
            if op not in ['refversion', 'reference']:
                print '%s: %s' % (op, config.get(rv, op))

    return


def add_parser_config(subparsers):

    parser = subparsers.add_parser('config', help=__doc__)
    parser.add_argument('-k', default=None, help='key')
    parser.add_argument('-v', default=None, help='value')
    parser.add_argument('-s', default='DEFAULT', help='reference version')
    parser.add_argument('--download_hg19', action='store_true', help='download hg19 reference and annotations')
    parser.add_argument('--download_hg18_anno', action='store_true', help='download hg18 (GRCh36) annotations')
    parser.add_argument('--download_hg19_anno', action='store_true', help='download hg19 (GRCh37) annotations')
    parser.add_argument('--download_hg38_anno', action='store_true', help='download hg38 (GRCh38) annotations')
    parser.add_argument('--download_mm10_anno', action='store_true', help='download mm10 (GRCm38) annotations')
    parser.add_argument('--download_mm9_anno', action='store_true', help='download mm9 (NCBIM37) annotations')
    parser.add_argument('--download_hg19_dbsnp', action='store_true', help='download hg19 dbsnp')
    parser.add_argument('--download_idmap', action='store_true', help='download id map')
    parser.set_defaults(func=main)

def add_parser_current(subparsers):

    parser = subparsers.add_parser('current', help="view current config")
    parser.set_defaults(func=main_current)
