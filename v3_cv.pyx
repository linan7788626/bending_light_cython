#cython: language_level=3,boundscheck=False, wraparound=False,nonecheck=False
#import cython
import numpy as np
cimport numpy as np

ctypedef np.float64_t Dtype
#ctypedef double Dtype

#-------------------------------------------------------------

cdef extern from "all_cv_test.h":
    void lens_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_lens)
    void mmbr_images(double *xi1,double *xi2,int nx1,int nx2,double *gpar,int npars,double *gpars,int nsubs,double *g_edge)
    void all_about_lensing(double *xi1,double *xi2,int nx1,int nx2,double * spar,double * lpar,int npars,double * lpars,int nsubs,double *s_image,double *g_lensimage,double *critical,double *caustic)

def call_lens_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                     np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                     np.ndarray[Dtype, ndim=1, mode="c"] gpar,
                     list gpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef int npars = len(gpar)
    cdef int nsubs = len(gpars)
    cdef np.ndarray gpars_array = np.zeros((nsubs,npars),dtype=np.double)
    gpars_array = np.array(gpars)
    cdef np.ndarray g_lens = np.zeros((nx1,nx2),dtype=np.double)
    lens_images(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpar.data,npars,<Dtype *>gpars_array.data,nsubs,<Dtype *>g_lens.data)

    return g_lens

def call_mmbr_images(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                     np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                     np.ndarray[Dtype, ndim=1, mode="c"] gpar,
                     list gpars):
    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef int npars = len(gpar)
    cdef int nsubs = len(gpars)
    cdef np.ndarray gpars_array = np.zeros((nsubs,npars),dtype=np.double)
    gpars_array = np.array(gpars)
    cdef np.ndarray g_edge = np.zeros((nx1,nx2),dtype=np.double)
    mmbr_images(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>gpar.data,npars,<Dtype *>gpars_array.data,nsubs,<Dtype *>g_edge.data)

    return g_edge

def call_all_about_lensing(np.ndarray[Dtype, ndim=2, mode="c"] xi1,
                      np.ndarray[Dtype, ndim=2, mode="c"] xi2,
                      np.ndarray[Dtype, ndim=1, mode="c"] spar,
                      np.ndarray[Dtype, ndim=1, mode="c"] lpar,
                      list lpars):

    cdef int nx1
    cdef int nx2
    nx1,nx2 = np.shape(xi1)
    cdef int npars = len(lpar)
    cdef int nsubs = len(lpars)
    cdef np.ndarray lpars_array = np.zeros((nsubs,npars),dtype=np.double)
    lpars_array = np.array(lpars)

    cdef np.ndarray s_image = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray g_lensimage = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray critical = np.zeros((nx1,nx2),dtype=np.double)
    cdef np.ndarray caustic = np.zeros((nx1,nx2),dtype=np.double)

    all_about_lensing(<Dtype *>xi1.data,<Dtype *>xi2.data,nx1,nx2,<Dtype *>spar.data,<Dtype *>lpar.data,npars,<Dtype *>lpars_array.data,nsubs,<Dtype *>s_image.data,<Dtype *>g_lensimage.data,<Dtype *>critical.data,<Dtype *>caustic.data)

    return s_image,g_lensimage,critical,caustic.T
