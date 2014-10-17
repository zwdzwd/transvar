""" configure TransVar """
import ConfigParser
import os

cfg_fns = [os.path.join(os.path.dirname(__file__), 'transvar.cfg'),
           os.path.expanduser('~/.transvar.cfg')]

downloaddirs = [os.path.join(os.path.dirname(__file__), 'transvar.download'),
                os.path.expanduser('~/.transvar.download')]

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
        if not os.path.exists(pdir): os.makedirs(pdir)

        try:
            for k, fn, url in fns:
                fnn = os.path.join(pdir, fn)
                download_url(url, fnn)
                if k:
                    config_set(config, section, k, fnn)
        except:
            continue

        break

def download_hg18_annotations(config):

    fns = [
        ('refseq', 'hg18.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/ARCHIVE/BUILD.36.3/GFF/ref_NCBI36_top_level.gff3.gz'),
        ('ccds', 'hg18.ccds.txt', 'ftp://ftp.ncbi.nlm.nih.gov/pub/CCDS/archive/Hs36.3/CCDS.20090327.txt'),
        ('aceview', 'hg18.aceview.gff.gz', 'ftp://ftp.ncbi.nih.gov/repository/acedb/ncbi_36_Apr07.human.genes/AceView.ncbi_36.genes_gff.tar.gz'),
        ('gencode', 'hg18.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_3c/gencode.v3c.annotation.NCBI36.gtf.gz'),
        ('ucsc', 'hg18.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg18.ucsc.refgene.txt.gz?dl=1'),
        ('ensembl', 'hg18.ensembl', 'ftp://ftp.ensembl.org/pub/release-54/gtf/homo_sapiens/Homo_sapiens.NCBI36.54.gtf.gz'),
        ]

    config.set('DEFAULT', 'refversion', 'hg18')
    _download_(config, 'hg18', fns)

def download_hg19_annotations(config):

    fns = [
        ('refseq', 'hg19.refseq.gff.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/ref_GRCh37.p13_top_level.gff3.gz?dl=1'),
        ('ccds', 'hg19.ccds.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/CCDS.current.txt?dl=1'),
        ('ensembl', 'hg19.ensembl.gtf.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/Homo_sapiens.GRCh37.75.gtf.gz?dl=1'),
        ('gencode', 'hg19.gencode.gtf.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/gencode.v19.annotation.gtf.gz?dl=1'),
        ('ucsc', 'hg19.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.ucsc.refgene.txt.gz?dl=1'),
        ('custom', 'hg19.custom.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.map?dl=1'),
        ('aceview', 'hg19.aceview.gff.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/AceView.ncbi_37.genes_gff.gff.gz?dl=1'),
        ('known_gene', 'hg19.knowngene.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_knownGene_hg19.gz?dl=1'),
        ('known_gene_alias', 'hg19.knowgene_alias.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_kgAlias.gz?dl=1'),
        ]

    config.set('DEFAULT', 'refversion', 'hg19')
    _download_(config, 'hg19', fns)

def download_idmap(config):
    fns = [('uniprot', 'uniprot.idmapping.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/HUMAN_9606_idmapping.dat.gz?dl=1')]
    _download_(config, 'idmap', fns)

def download_hg38_annotations(config):

    fns = [('refseq', 'hg38.refseq.gff.gz', 'ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38_top_level.gff3.gz'),
           # ('ccds', 'hg38.ccds.txt', ''),
           ('ensembl', 'hg38.ensembl.gtf.gz', 'ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz'),
           ('gencode', 'hg38.gencode.gtf.gz', 'ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz'),
           ('ucsc', 'hg38.ucsc.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg38.ucsc.refgene.txt.gz?dl=1'),
       ]
    config.set('DEFAULT', 'refversion', 'hg38')
    _download_(config, 'hg38', fns)

# def download_hg19_reference(config):
#     fns = [('reference', 'hg19.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)
    
# def download_hg38_reference(config):
#     fns = [('reference', 'hg38.fa.gz', 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz')]
#     _download_(config, 'hg38', fns)

def download_hg19_dbsnp(config):
    fns = [('dbsnp', 'hg19_dbsnp.vcf.gz', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz'),
           ('dbsnp_index', 'hg19_dbsnp.vcf.gz.tbi', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz.tbi'),
    ]

    _download_(config, 'hg19', fns)

def read_config():
    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    return config

def main(args):

    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    if args.k and args.v:
        config_set(config, args.s, args.k, args.v)

    if args.download_hg18_anno:
        download_hg18_annotations(config)

    if args.download_hg19_anno:
        download_hg19_annotations(config)

    if args.download_hg38_anno:
        download_hg38_annotations(config)

    if args.download_hg19_dbsnp:
        download_hg19_dbsnp(config)

    if args.download_idmap:
        download_idmap(config)

    for cfg_fn in cfg_fns:
        try:
            config.write(open(cfg_fn,'w'))
            break
        except IOError as e:
            pass


def add_parser_config(subparsers):

    parser = subparsers.add_parser('config', help=__doc__)
    parser.add_argument('-k', default=None, help='key')
    parser.add_argument('-v', default=None, help='value')
    parser.add_argument('-s', default='DEFAULT', help='reference version')
    parser.add_argument('--download_hg19', action='store_true', help='download hg19 reference and annotations')
    parser.add_argument('--download_hg18_anno', action='store_true', help='download hg18 (GRCh36) annotations')
    parser.add_argument('--download_hg19_anno', action='store_true', help='download hg19 (GRCh37) annotations')
    parser.add_argument('--download_hg38_anno', action='store_true', help='download hg38 (GRCh38) annotations')
    parser.add_argument('--download_hg19_dbsnp', action='store_true', help='download hg19 dbsnp')
    parser.add_argument('--download_idmap', action='store_true', help='download id map')
    parser.set_defaults(func=main)
