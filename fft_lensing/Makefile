#for Gnu's gcc compiler, -O2 for optimization, -g for debugging
CC        = gcc
CC_OPTS_GSL = -lgsl -lgslcblas  
CC_OPTS_FFTW  =  -lfftw3
CC_OPTS_CFITS =  -lcfitsio 
#CC_OPTS_LENSING = -I../include -L../lib -lgglensing 
CC_OPTS_LENSING = -I../include -L../lib -lglsg 

CC_WARNINGS  = -O2 -Wall -lm

all: cal_alphas

clean: 
	rm cal_alphas 

cal_alphas: main.c
	$(CC) main.c $(CC_OPTS_LENSING) $(CC_OPTS_GSL) $(CC_OPTS_FFTW) $(CC_OPTS_CFITS) $(CC_WARNINGS) -o cal_alphas 

# end of Makefile
