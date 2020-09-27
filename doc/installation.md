#Installation
##Before you start
SPLICE-q requires Python 3.6+. To check your Python version, run in your terminal (Mac/Linux/Win):
```bash
 $ python --version
```
You should get a output like `Python 3.6.3`. If you do not have Python (or need to update it), please install the latest 3.x version from [python.org](https://www.python.org/downloads/).

## Use pip for installing
Using pip is the easiest way to install SPLICE-q. First, make sure you have [pip](https://packaging.python.org/key_projects/#pip) available.

```bash
 $ pip --version
```

You should see an output displaying the pip version, as well as the location and Python version. If you don't, please install it by following the instructions for your system as described [here](https://pip.pypa.io/en/stable/installing/).

##Install SPLICE-q
SPLICE-q can be installed from pip and from source.
### pip

```bash
 $ pip3 install SPLICE-q
```

###Development/install from source

```bash
 $ git clone https://github.com/vrmelo/SPLICE-q
 $ cd SPLICE-q
 $ pip3 install -e .
```

Requirements

- Python 3.6+
- PySam
- InterLap
- NumPy
- Rich

Operating Systems
- Linux, macOS, and Windows 10 Subsystem for Linux.
