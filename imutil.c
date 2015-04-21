/*
 * Author: Francois Rigaut
 *
 * This file contains a number of utility functions, coded in C to gain
 * execution time. It addresses functionalities that are missing in
 * yorick, mostly concerning 2D image processing.
 *
 * Copyright (c) 2003-2011, Francois Rigaut
 *
 * This program is free software; you can redistribute it and/or  modify it
 * under the terms of the GNU General Public License  as  published  by the
 * Free Software Foundation; either version 2 of the License,  or  (at your
 * option) any later version.
 *
 * This program is distributed in the hope  that  it  will  be  useful, but
 * WITHOUT  ANY   WARRANTY;   without   even   the   implied   warranty  of
 * MERCHANTABILITY or  FITNESS  FOR  A  PARTICULAR  PURPOSE.   See  the GNU
 * General Public License for more details (to receive a  copy  of  the GNU
 * General Public License, write to the Free Software Foundation, Inc., 675
 * Mass Ave, Cambridge, MA 02139, USA).
 *
 */


#include <math.h>
#include <stdio.h>

void srandom(unsigned seed);


/************************************************************************
 * noop. For testing and timing.                                        *
 ************************************************************************/

int _mynoop1()
{
  return (0);
}

void _dpc(float *ima, short *gpm, long *w, long *w2, long stride, long npix, long nw, long npassmax)
{
  long i,k,n,npass;
  float val;

  npass = 0;

  while (nw && (npass<npassmax)) {

    for (i=0;i<nw;i++) {
      val = 0.0f;
      n = 0l;
      w2[i] = 0l;
      k = w[i]-1;      if ((k>=0)&&(k<npix)&&(gpm[k])) { val += ima[k]; n++; }
      k = w[i]+1;      if ((k>=0)&&(k<npix)&&(gpm[k])) { val += ima[k]; n++; }
      k = w[i]-stride; if ((k>=0)&&(k<npix)&&(gpm[k])) { val += ima[k]; n++; }
      k = w[i]+stride; if ((k>=0)&&(k<npix)&&(gpm[k])) { val += ima[k]; n++; }
      if (n) {
        // fill in corrected value in image:
        ima[w[i]] = val/(float)n;
        // mark the pixel as processed:
        w2[i] = 1;
      }
    }
    k = 0;
    for (i=0;i<nw;i++) {
      if (w2[i]) { // pixel has been processed:
        gpm[w[i]] = 1;
      } else { // pixel is still bad:
        w[k++] = w[i];
      }
    }
    nw = k;
    npass++;
  }
}

void _sigfil(float *ima, short *bpm, float nsigma, long stride, long npix)
{
  long i,k,n;
  float avg1,avg2,val,std;

  for (i=0;i<npix;i++) {
    avg1 = avg2 = 0.0f;
    n = 0;
    // compute the avg and rms in a 3x3 box centered on pixel
    // k = i; avg1=ima[k]; avg2=ima[k]*ima[k]; n++;
    k = i-1;        if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+1;        if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride;   if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride-1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i-stride+1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride;   if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride-1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    k = i+stride+1; if ((k>=0)&&(k<npix)) { avg1+=ima[k]; avg2+=ima[k]*ima[k]; n++; }
    avg1 = avg1 / (float)n;
    avg2 = avg2 / (float)n;
    std  = avg2 - avg1 * avg1;
    if (std<0.) std=0.0f;
    else std = sqrtf(std);
    val = fabsf(ima[i]-avg1);
    // printf("ima=%f val=%f avg1=%f, avg2=%f, std=%f nsigma=%f\n",ima[i],val,avg1,avg2,std,nsigma);
    if (val>(nsigma*std)) bpm[i]=1; else bpm[i]=0;
  }
}

/************************************************************************
 * Functions _bin2d                                                     *
 * Returns the input image, rebinned with the specified binning factor  *
 * The input image dimension needs not to be a multiple of the binning  *
 * factor. If it is not the case, the final pixel counts several times  *
 * the edge pixels.                                                     *
 * There are 3 function varieties to deal with longs, float and double. *
 * Called by function bin2d(in,binfact) in yorickUtils.i                *
 * Last modified: December 15, 2003.                                    *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

int _bin2d_long(long *in, int nx, int ny, long *out, int fx, int fy, int binfact)
{
  /* Declarations */
  int i1,i2,j1,j2,i,j;

  /* Loop on indices to bin */
  for ( i1=0 ; i1<fx ; i1++ ) {
    for ( j1=0 ; j1<fy ; j1++ ) {
      for ( i2=0 ; i2<binfact ; i2++ ) {
  for ( j2=0 ; j2<binfact ; j2++ ) {
    i = i1*binfact+i2; if ( i>=nx ) { i=nx-1;}
    j = j1*binfact+j2; if ( j>=ny ) { j=ny-1;}
    out[i1+j1*fx] += in[i+j*nx];
  }
      }
    }
  }
  return (0);
}

int _bin2d_float(float *in, int nx, int ny, float *out, int fx, int fy, int binfact)
{
  /* Declarations */
  int i1,i2,j1,j2,i,j;

  /* Loop on indices to bin */
  for ( i1=0 ; i1<fx ; i1++ ) {
    for ( j1=0 ; j1<fy ; j1++ ) {
      for ( i2=0 ; i2<binfact ; i2++ ) {
  for ( j2=0 ; j2<binfact ; j2++ ) {
    i = i1*binfact+i2; if ( i>=nx ) { i=nx-1;}
    j = j1*binfact+j2; if ( j>=ny ) { j=ny-1;}
    out[i1+j1*fx] += in[i+j*nx];
  }
      }
    }
  }
  return (0);
}

int _bin2d_double(double *in, int nx, int ny, double *out, int fx, int fy, int binfact)
{
  /* Declarations */
  int i1,i2,j1,j2,i,j;

  /* Loop on indices to bin */
  for ( i1=0 ; i1<fx ; i1++ ) {
    for ( j1=0 ; j1<fy ; j1++ ) {
      for ( i2=0 ; i2<binfact ; i2++ ) {
  for ( j2=0 ; j2<binfact ; j2++ ) {
    i = i1*binfact+i2; if ( i>=nx ) { i=nx-1;}
    j = j1*binfact+j2; if ( j>=ny ) { j=ny-1;}
    out[i1+j1*fx] += in[i+j*nx];
  }
      }
    }
  }
  return (0);
}


/************************************************************************
 * Function _eclat                                                      *
 * Returns results identical to roll(), but faster, as it is dedicated  *
 * to swapping quadrants, for use with FFTs.                            *
 * Warning: In-place swapping.
 * Last modified: December 15, 2003.                                    *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

void _eclat_long(long *ar, int nx, int ny)
{
  int i,j,k1,k2;
  long a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}

void _eclat_float(float *ar, int nx, int ny)
{
  int i,j,k1,k2;
  float a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}

void _eclat_double(double *ar, int nx, int ny)
{
  int i,j,k1,k2;
  double a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}

/************************************************************************
 * Function _dist                                                       *
 * returns a 2D array which element's values contain the distance to    *
 * the point of coordinates (xc,yc).                                    *
 * Called by dist() in yorickfr.i                                       *
 * Last modified: December 15, 2003.                                    *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

void _dist(float *d, long dimx, long dimy, float xc, float yc)
{
  /* Declarations */
  long i,j;

  /* Loop and fill d with distance values */
  for (i=0;i<dimx;++i) {
    for (j=0;j<dimy;++j) {
      d[i + j * dimx] = (float)sqrt( (xc-(float)i) * (xc-(float)i) +
             (yc-(float)j) * (yc-(float)j) );
    }
  }
}


/************************************************************************
 * Functions _clip                                                      *
 * Returns the input array, in which elements smaller than xmin have    *
 * been replaced by xmin, and greater than xmax by xmax.                *
 * There is a serie of these routines to handle long, float and double  *
 * types. The clipmin and clipmax is to clip only the min or only the   *
 * max value.
 * Called by clip() in yorickfr.i                                       *
 * Last modified: December 15, 2003.                                    *
 * Author: F.Rigaut                                                     *
 ************************************************************************/

int clipchar(char *x, char xmin, char xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipshort(short *x, short xmin, short xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipint(int *x, int xmin, int xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int cliplong(long *x, long xmin, long xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipfloat(float *x, float xmin, float xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipdouble(double *x, double xmin, double xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

/***********************/

int clipminchar(char *x, char xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

int clipminshort(short *x, short xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

int clipminint(int *x, int xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

int clipminlong(long *x, long xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

int clipminfloat(float *x, float xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

int clipmindouble(double *x, double xmin, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] < xmin) x[i]=xmin;
  }
  return (0);
}

/*********************/

int clipmaxchar(char *x, char xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipmaxshort(short *x, short xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipmaxint(int *x, int xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipmaxlong(long *x, long xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipmaxfloat(float *x, float xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

int clipmaxdouble(double *x, double xmax, long n)
{
  long i = 0;
  for (i=0;i<n;++i) {
    if (x[i] > xmax) x[i]=xmax;
  }
  return (0);
}

void _ran1init(long n)
{
  srandom(n);
}

float ran1()
{
  float norm;
  norm = 2147483647.f;

  return random()/norm;
}

void _poidev(float *xmv, long n)
/* all floats -> doubles on June 2010 to avoid SIGFPE
   for too large input values */
{
  double gammln(double xx);
  /*  float ran1(long *idum);*/
  static double sq,alxm,g,oldm=(-1.0);
  double xm,em,t,y,y1;
  long i;

  for (i=0;i<n;i++) {
    xm = (double)xmv[i];
    if (xm == 0.0f) continue;
    if (xm < 20.0) { /* Use direct method. */
      if (xm != oldm) {
        oldm=xm;
        g=exp(-xm);  /* If xm is new, compute the exponential. */
      }
      em = -1;
      t=1.0;
      do {
        ++em;
        t *= ran1();
      } while (t > g);
    } else {  /* Use rejection method. */
      if (xm != oldm) {
        oldm=xm;
        sq=sqrt(2.0*xm);
        alxm=log(xm);
        // printf("xm+1.0 = %.f gammln(xm+1.0) = %.f\n",xm+1.0,gammln(xm+1.0));
        g=xm*alxm-gammln(xm+1.0);
      }
      do {
        do {
          y=tan(3.1415926535897932384626433832*ran1());
          em=sq*y+xm;
        } while (em < 0.0);
        em=floor(em);
        // printf("em+1.0 = %.f gammln(em+1.0) = %.f\n",em+1.0,gammln(em+1.0));
        // printf("exp(em*alxm-gammln(em+1.0)-g) = %.f\n",exp(em*alxm-gammln(em+1.0)-g));
        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (ran1() > t);
    }
    xmv[i] = (float)em;
  }
}


double gammln(double xx)
{
  /* Returns the value ln[?(xx)] for xx>0. */
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
      24.01409824083091,-1.231739572450155,
      0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}


void _gaussdev(float *xmv, long n)
{
  /* Returns a normally distributed deviate with zero mean and unit variance,
     using ran1() as the source of uniform deviates. */

  /*  float ran1(long *idum); */
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  long i;

  for (i=0;i<n;i++) {
    if (iset == 0) {
      do {
  v1=2.0*ran1()-1.0;
  v2=2.0*ran1()-1.0;
  rsq=v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      gset=v1*fac;
      iset=1;
      xmv[i] = v2*fac;
    } else {
      iset=0;
      xmv[i] = gset;
    }
  }
}

