
# TransVar [![Travis-CI Build Status](https://travis-ci.org/zwdzwd/transvar.svg?branch=master)](https://travis-ci.org/zwdzwd/transvar) [![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org) [![PyPI Downloads](https://img.shields.io/pypi/dm/transvar.svg)](https://pypi.org/project/TransVar/)

TransVar is a multi-way annotator for genetic elements and genetic variations. It operates on genomic coordinates (e.g., `chr3:g.178936091G>A`) and transcript-dependent cDNA as well as protein coordinates (e.g., `PIK3CA:p.E545K` or `PIK3CA:c.1633G>A`, or `NM_006218.2:p.E545K`, or `NP_006266.2:p.G240Afs*50`). It is particularly designed with the functionality of resolving ambiguous mutation annotations arising from differential transcript usage. TransVar keeps awareness of the underlying unknown transcript structure (exon boundary, reference amino acid/base) while performing reverse annotation (via fuzzy matching from protein level to cDNA level).

**User Guide is [here](http://transvar.readthedocs.io/en/latest/).**

Try out transvar via pip

```bash
sudo pip install transvar
```
or locally
```bash
pip install --user transvar
```

To upgrade from previous installation, you can
```bash
pip install -U transvar
```

Pre-built docker images can be found [here](https://cloud.docker.com/repository/docker/zhouwanding/transvar/general)

This is a continued TransVar implementation from what was hosted at [https://bitbucket.org/wanding/transvar](https://bitbucket.org/wanding/transvar).
