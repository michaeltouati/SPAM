[<img src='https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54' height="20">](https://www.python.org/)
[<img src='https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white' height="20">](https://numpy.org/)
[<img src='https://matplotlib.org/_static/logo2_compressed.svg' height="20">](https://matplotlib.org/stable/index.html#)
[<img src='https://mpng.subpng.com/20180408/xze/kisspng-wing-ide-integrated-development-environment-python-raspberry-5aca9cd85768a8.1913876415232278643581.jpg' height="20">](https://docs.python.org/3/library/tkinter.html)
[<img src='https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white' height="20">](https://www.latex-project.org//)
[![License : GPLv3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)
[![Downloads](https://img.shields.io/github/downloads/michaeltouati/SPAM/total)](https://github.com/michaeltouati/SPAM/releases)
[![Views : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Views%20(<15%20days)&query=count&url=https://github.com/michaeltouati/SPAM/blob/master/.github/view.json?raw=True&logo=github)](https://github.com/michaeltouati/SPAM/actions/workflows/views.yml)
[![Clones : 14 days](https://img.shields.io/badge/dynamic/json?color=success&label=Clones%20(<15%20days)&query=count&url=https://github.com/michaeltouati/SPAM/blob/master/.github/clone.json?raw=True&logo=github)](https://github.com/michaeltouati/SPAM/actions/workflows/clones.yml)
[![Open issues](https://img.shields.io/github/issues/michaeltouati/SPAM)](https://github.com/michaeltouati/SPAM/issues)
[![Closed issues](https://img.shields.io/github/issues-closed/michaeltouati/SPAM)](https://github.com/michaeltouati/SPAM/issues)
[![Open pull requests](https://img.shields.io/github/issues-pr/michaeltouati/SPAM)](https://github.com/michaeltouati/SPAM/pulls)
[![Closed pull requests](https://img.shields.io/github/issues-pr-closed/michaeltouati/SPAM)](https://github.com/michaeltouati/SPAM/pulls)

# SPAM
Written by Michaël TOUATI  - 09/27/2019

<p align="center">
  <img width="800" height="450" src="spam.png">
</p>

SPAM (Stopping Power of Protons and Alpha particles in Ambiant Matter) is a python tool allowing for printing and saving the stopping power and/or the Bragg's peak of protons or alpha particles in ambiant matter. Equations computed by SPAM are detailed in the spam.pdf file. The comparisons with [PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) and [ASTAR](https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html) NIST databases that has been used to validate SPAM can be found in the benchmarks.tar.xz file. 

SPAM needs the following Python packages :
- tkinter
- PIL 
- numpy 
- matplotlib

The next release will take into account the angular scattering of protons and alpha particles using a Monte-Carlo approach. 

In order to use it, download or clone it and just type

```sh
python3 spam.py
```

The tkinter window will pop up.

