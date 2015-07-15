# cython_pylensing_toys

=======================================
## Introduction
A version of my glensing toy with cython boosting.

###Dependences
Python2.7, Numpy, pygame, cython.

###Compile
"python setup.py build_ext --inplace" OR "./make_so"

###Run test
python test_cython.py


##Usage

###How to control sources.
Use the left click of the mouse to control properties of source, 
but you have to press "w" for where, "s" for size and "e" for ellipticity and orientation.
* Holding "w" and left click, you can move you mouse to move the position of the source.
* Holding "s" and left click, you can move you mouse up and down to change the size of the source.
* Holding "e" and left click, you can move you mouse up and down to change the ellipticity of the source, or, left and right to change the orientation of the source.

###how to control lenses: 
Use the right click of the mouse to control properties of lenses, 
as controling the sources, you can also press "w" for position, "s" for lensing strength and "e" for ellipticity and orientation, but first of all, you need to press num keys to choose which component you want to handle. Here, "1" stands for the main halo, "2~9" stand for the sub halos. For example, if you want to change the ellipticity and orientation of the main halo, please hold "1", "e" and right click, then move vertically, the ellipticity of the main halo will be changed, move horizontally, the orientation of the main halo will be changed. Because we are using NIE model in this code, so when you press "s" to change the size of the lens, there are to sizes, one the Einstein Radius(move your mouse vertically), the other one is the core radius (move your mouse horizontally). If you want to control the subhalos, for example, at the beginning, there is only one subhalo, please press "2" first, then do everythin you want.

###How to add or remove satellite galaxies.
In our code, we can also add or remove subhalos as you will, but, so far, you can only add 8 subhalos at most.
>Hold "=", click the right click to add.
>Hold "-", click the right click to remove the last added subhalo.

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
