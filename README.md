# SPAM
Written by Dr M Touati - CLPU - September 2019 - 
mtouati@clpu.es

SPAM (Stopping Power of Protons and Alpha particles in Ambiant Matter) is a 
python tool allowing for printing and saving the stopping power and/or the 
Bragg's peak of protons or alpha particles in ambiant matter.
You can find the litterature that has been used to write the code in the 
folder litterature as well as the benchmarks performed by comparing SPAM with NIST databases (https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html and https://physics.nist.gov/PhysRefData/Star/Text/ASTAR.html) in the tar file benchmarks.tar.xz. 

It needs the following Python packages :
- tkinter
- PIL 
- numpy 
- matplotlib

The next release will take into account the angular scattering of protons and alpha particles using a Monte-Carlo approach. 

Also, in order to avoid the code to be slow because of the search by dichotomy of values in tabulated functions and tables taken from the litterature, it is planned to use Cython instead.

In order to use it download it and just type

$python3 spam.py

A Tkinter window will open and you can chose the projectile as well as the material.

