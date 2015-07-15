# cython_pylensing_toys

=======================================
## Introduction
A version of my glensing toy with cython boosting. 

1. Dependences: 
>Python2.7, Numpy, pygame, cython.

2. Installation:
>* Clone it to your local directory or click the "Download ZIP" button to download it.  
>* "python setup.py build_ext --inplace" OR "./make_so"

3. Run Code:
>python test_cython.py


##Usage

1. How to control source:  
Use the left click of the mouse to control the source, 
but you have to press "w" for where, "s" for size and "e" for ellipticity and orientation.
  * Holding "w" and left click, you can move you mouse to move the position of the source.
  * Holding "s" and left click, you can move you mouse up and down to change the size of the source.
  * Holding "e" and left click, you can move you mouse up and down to change the ellipticity of the source, or, left and right to change the orientation of the source.

2. How to control lenses:  
Use the right click of the mouse to control the lenses.
Similar with the source, you can press "w" for position, "s" for size (Einstein radius and core radius) and "e" for ellipticity and orientation. But first of all, you need to press one num key to choose a component you want to handle. Here, "1" stands for the main halo, "2~9" stand for the sub halos. For example, if you want to change the ellipticity and orientation of the main halo: please hold "1", "e" and right click, then move vertically, the ellipticity of the main halo will be changed, move horizontally, the orientation of the main halo will be changed. Because we are using NIE model in this code, so when you press "s" to change the size of the lens, there are two sizes, one the Einstein Radius(move your mouse vertically), the other one is the core radius (move your mouse horizontally). If you want to control the subhalos, there is only one subhalo intrinsicly, please press "2" first, then do everything you want.

  * Holding "w" and right click, you can move you mouse to move the position of the chosen component of the lens.
  * Holding "s" and right click, you can move you mouse vertically to change the Einstein radius of chosen component, or, horizontally to change the Core radius of the chosen component.
  * Holding "e" and right click, you can move you mouse vertically  to change the ellipticity of the chosen component, or, horizontally to change the orientation of the chosen component.

3. How to add or remove satellite galaxies:  
In our code, we can also add or remove subhalos if needed, but, so far, you can only add 8 subhalos at most.
  *Hold "=", click the right click to add.
  *Hold "-", click the right click to remove the last added subhalo.

4. Toggles:  
Try "f" and "g". :-)

##Todo List
1. Sources
 - [ ] More models of sources (Sersic, Moffat, Disk\+Bulge...)
 - [ ] Real galaxies images.
 - [ ] Arbitrary pictures (for fun).

2. Lenses
 - [ ] More models of lenses (PIEMD, GNFW, NFW\+Hernquist, Burkert...)
 - [ ] Directly input mass sheet.
 - [ ] More satellite galaxies?
 - [ ] Realistic images of statelite galaxies.

3. Algorithm 
 - [ ] Better caustic drawer.

5. Time delay
 - [ ] Implementing in twinkles branch.

6. Competition Model.
 - [ ] Input lensed images by one person.
 - [ ] Reconstruct it by another one.
 - [ ] Score the matching.

7. Science Side
 - [ ] Save the parameters of current lensing system, including lenses and sources.
 - [ ] Make the catalog of lensing system to be an input file.
 - [ ] Show all the parameters of current lensing system in real time.
