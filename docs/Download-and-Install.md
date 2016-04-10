## Install using pip

```bash
sudo pip install transvar
```
or locally
```bash
pip install --user transvar
```

to upgrade from a previous version
```bash
pip install -U transvar
```

## Download the program

Latest release is available [here](https://github.com/zwdzwd/transvar/releases/latest)

For all previous versions, see [here](https://github.com/zwdzwd/transvar/releases).

### Other old stable releases
 + stable 2.0.x version [v2.0.12.20150626](https://github.com/zwdzwd/transvar/archive/v2.0.12.20150626.zip)
 + stable 1.x version [v1.40](https://github.com/zwdzwd/transvar/archive/v1.40.zip)

## Dependency

The only requirement for building TransVar are Python 2.7 and a reasonably modern C compiler such as gcc.

## Install from source

### Local install
```bash
python setup.py install --prefix [folder]
```
The installation will create two subfolders of `[folder]`: `[folder]/lib` (which would contain libraries) and `[folder]/bin` (which would contain transvar executable).
When you run transvar, make sure `[folder]/lib/python2.7/site-packages` is in your PYTHONPATH. In some occasions, you need to `mkdir -p [folder]/lib/python2.7/site-packages` to make sure it exists before you could run `setup.py`.
You can add it by putting
`export PYTHONPATH=$PYTHONPATH:[folder]/lib/python-2.7/site-packages/` to your `.bashrc` (or `.profile` depending on your OS).

The installed executable is `[folder]/bin/transvar`.

### System-wise install (need root)
```bash
sudo python setup.py install
```
