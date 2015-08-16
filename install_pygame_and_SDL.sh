# make a temporary directory where we can download and build stuff
mkdir tmp
cd tmp

# download and install SDL
wget http://www.libsdl.org/release/SDL-1.2.14.tar.gz
tar -xzvf SDL-1.2.14.tar.gz
cd SDL-1.2.14
./configure --prefix=$HOME
make
make install

# download and extract PyGame
wget http://pygame.seul.org/ftp/pygame-1.9.1release.tar.gz
tar xzvf pygame-1.9.1release.tar.gz
cd pygame-1.9.1release

# Here you need to edit the Setup file and comment out the line that looks like
# _camera src/_camera.c src/camera_v4l2.c src/camera_v4l.c $(SDL) $(DEBUG)
# for some reason the compilation of the camera module fails for me, so I commented it out
vi Setup

# once the camera module is out of the way we can proceed with the installation
# Make sure to use the correct Python version here (e.g. 'python2.5' or 'python2.6')
python2.5 setup.py install --prefix=$HOME

# That's it


####################################################
To get a clean compile, edit ./src/video/x11/SDL_x11sym.h around line 166 so it looks like this:

Code:

#if 0
#ifdef LONG64
SDL_X11_MODULE(IO_32BIT)
SDL_X11_SYM(int,_XData32,(Display *dpy,register long *data,unsigned len),(dpy,data,len),return)
SDL_X11_SYM(void,_XRead32,(Display *dpy,register long *data,long len),(dpy,data,len),)
#endif
#endif
####################################################

