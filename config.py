""" configure revan """
import ConfigParser
import os

cfg_fns = [os.path.join(os.path.dirname(__file__), 'revan.cfg'),
           os.path.expanduser('~/.revan.cfg')]

def download_hg19_annotations(config):

    pass

def download_hg19_reference(config):

    pass

def read_config():
    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    return config.defaults()    

def main(args):

    config = ConfigParser.RawConfigParser()
    config.read(cfg_fns)
    config.set('DEFAULT', args.k, args.v)
    d = config.defaults()

    if args.download_hg19:
        download_hg19_annotations(config)
        download_hg19_reference(config)

    if args.download_hg19_anno:
        download_hg19_annotations(config)

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
    parser.set_defaults(func=main)
