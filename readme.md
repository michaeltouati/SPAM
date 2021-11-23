# SPAM
Written by Michaël TOUATI  - 09/27/2019

SPAM (Stopping Power of Protons and Alpha particles in Ambiant Matter) is a python tool allowing for printing and saving the stopping power and/or the Bragg's peak of protons or alpha particles in ambiant matter. Equations computed by SPAM are detailed in the spam.pdf file. All benchmarks performed by comparing SPAM with [PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) and [ASTAR](https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html) NIST databases can be found in the benchmarks.tar.xz file . 

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

