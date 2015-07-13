# cython_pylensing_toys

=======================================
## Introduction
A version of my glensing toy with cython boosting.

###Dependences
Python2.7, Numpy, pygame, cython.

###Compile
python setup.py build_ext --inplace

###Run test
python test_cython.py


##Usage

###How to control sources.
Use the left click of the mouse to control properties of sources, 
but you have to press "w" for where, "s" for size and "e" for ellipticity and orientation.

###how to control lenses: 
Use the left click of the mouse to control properties of sources, 
but you have to press "w" for position, "s" for lensing strength and "e" for ellipticity and orientation.

###How to add or remove satellite galaxies.
Hit and hold "=", click the right click.
Hit and hold "-", click the right click.

###How to add or remove more sources.
Hit and hold "=", click the right click.
Hit and hold "-", click the right click.

###Toggles 
Try "f" and "g". :-)

##Todo List
###Sources
1. More models of sources (Sersic, Moffat, Disk\+Bulge...)
2. Real galaxies images.
3. Arbitrary pictures.

###Lenses
1. More models of lenses (PIEMD, GNFW, NFW\+Hernquist, Burkert...)
2. Directly input mass sheet.
3. More satellite galaxies?

###Algorithm 
1. Better caustic drawer.

###Time delay
1. Implement in twinkles branch.

###Competition Models.
1. Input lensed images by one person.
2. Reconstruct it by another one.
3. Score the matching.


###Science Side
1. Save the parameters of current lensing system, including lenses and sources.
2. Make the catalog of lensing system to be an input file.
3. Show all the parameters of current lensing system in real time.
