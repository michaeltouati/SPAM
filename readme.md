# SPAM
Written by Michaël TOUATI 

SPAM (Stopping Power of Protons and Alpha particles in Ambiant Matter) is a 
python tool allowing for printing and saving the stopping power and/or the 
Bragg's peak of protons or alpha particles in ambiant matter.
You can find the litterature that has been used to write the code in the 
folder litterature as well as the benchmarks performed by comparing SPAM with [PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html) and [ASTAR](https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html) NIST databases  in the tar file benchmarks.tar.xz. 

It needs the following Python packages :
- tkinter
- PIL 
- numpy 
- matplotlib

The next release will take into account the angular scattering of protons and alpha particles using a Monte-Carlo approach. 

In order to use it download or clone it and just type

$python3 spam.py

The tkinter window will pop up.

