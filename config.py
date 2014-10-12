""" configure TransVar """
import ConfigParser
import os

cfg_fns = [os.path.join(os.path.dirname(__file__), 'transvar.cfg'),
           os.path.expanduser('~/.transvar.cfg')]

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


def _download_(config, fns):

    pdir = os.path.join(os.path.dirname(__file__), 'download')
    if not os.path.exists(pdir): os.makedirs(pdir)

    for k, fn, url in fns:
        fnn = os.path.join(pdir, fn)
        download_url(url, fnn)
        config.set('DEFAULT', k, fnn)


def download_hg19_annotations(config):

    fns = [
        ('refseq', 'hg19.refseq.gff.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/ref_GRCh37.p13_top_level.gff3.gz?dl=1'),
        ('ccds', 'hg19.ccds.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/CCDS.current.txt?dl=1'),
        ('ensembl', 'hg19.ensembl.gtf.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/Homo_sapiens.GRCh37.75.gtf.gz?dl=1'),
        ('gencode', 'hg19.gencode.gtf.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/gencode.v19.annotation.gtf.gz?dl=1'),
        ('ucsc', 'hg19.ucsc.txt', 'https://dl.dropboxusercontent.com/u/6647241/annotations/hg19.map?dl=1'),
        ('aceview', 'hg19.aceview.gff.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/AceView.ncbi_37.genes_gff.gff.gz?dl=1'),
        ('known_gene', 'hg19.knowngene.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_knownGene_hg19.gz?dl=1'),
        ('known_gene_alias', 'hg19.knowgene_alias.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/UCSC_kgAlias.gz?dl=1'),
        ('uniprot', 'uniprot.idmapping.txt.gz', 'https://dl.dropboxusercontent.com/u/6647241/annotations/HUMAN_9606_idmapping.dat.gz?dl=1'),
        ]

    _download_(config, fns)

def download_hg19_reference(config):
    fns = ['reference', 'hg19_reference.fa', '']
    _download_(confg, fns)

def download_hg19_dbsnp(config):
    fns = [('dbsnp', 'hg19_dbsnp.vcf.gz', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz'),
           ('dbsnp_index', 'hg19_dbsnp.vcf.gz.tbi', 'ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b141_GRCh37p13/VCF/00-All.vcf.gz.tbi'),
    ]

    _download_(config, fns)

def read_config():
    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    return config.defaults()    

def main(args):

    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    if args.k and args.v:
        config.set('DEFAULT', args.k, args.v)
    d = config.defaults()

    if args.download_hg19:
        download_hg19_annotations(config)
        download_hg19_reference(config)

    if args.download_hg19_anno:
        download_hg19_annotations(config)

    if args.download_hg19_dbsnp:
        download_hg19_dbsnp(config)

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
    parser.add_argument('--download_hg19', action='store_true', help='download hg19 reference and annotations')
    parser.add_argument('--download_hg19_anno', action='store_true', help='download hg19 annotations')
    parser.add_argument('--download_hg19_dbsnp', action='store_true', help='download hg19 dbsnp')
    parser.set_defaults(func=main)
