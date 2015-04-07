def opengz(fn):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn)
    else:
        fh = open(fn)

    return fh
