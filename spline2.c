/*
 * Author: Francois Rigaut
 *
 * Copyright (c) 2003-2007, Francois Rigaut
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
 *
*/   

#include <stdio.h>
#include "ydata.h"
#include "pstdlib.h"

void _splint(float *xa, 
	     float *ya, 
	     float *y2a, 
	     long n, 
	     float x, 
	     float *y)
{
  long klo,khi,k;
  float h,b,a;
  
  klo=0;
  khi=n-1;
  while (khi-klo > 1) {
    k=(khi+klo) >> 1;
    if (xa[k] > x) khi=k;
    else klo=k;
  }
  h=xa[khi]-xa[klo];
  if (h == 0.0) YError("Bad xa input to routine _splint");
  a=(xa[khi]-x)/h;
  b=(x-xa[klo])/h;
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


void _splinf(float *x, 
	     float *y, 
	     long n, 
	     float *y2)
{
  long i,k;
  float p,qn,sig,un,*u;

  u = p_malloc(sizeof(float)*(n-1));
  
  y2[0]=u[0]=qn=un=0.0;
  
  for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  
  for (k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];
  
  p_free(u);
}

void 
_splin2(float *xin, 
	float *yin, 
	float *image, 
	float *deriv, 
	long nx, 
	long ny, 
	long *nvalidx, 
	float xout, 
	float yout, 
	float *res)
{
  long j,m;
  long n=0;
  float *y2tmp,*yytmp;

  y2tmp = p_malloc(sizeof(float)*ny);
  yytmp = p_malloc(sizeof(float)*ny);

  for (j=0;j<=ny-1;j++) {
    m = nvalidx[j];
    _splint(&xin[n],&image[n],&deriv[n],m,xout,&yytmp[j]);
    n += m;
  }

  _splinf(yin,yytmp,ny,y2tmp);
  _splint(yin,yytmp,y2tmp,ny,yout,res);

  p_free(y2tmp);
  p_free(yytmp);
}


void 
_splie2( float *xin, 
	 float *im, 
	 long nx, 
	 long ny,				\
	 float *deriv, 
	 long *nvalidx)
     /* updated for XY faster indice */
{
  long j,m;
  long n=0;

  for (j=0;j<=ny-1;j++) {
    m = nvalidx[j];
    _splinf(&xin[n],&im[n],m,&deriv[n]);
    n += m;
  }
}


void 
_spline2(  float *xin, 
	   float *yin, 
	   float *im, 
	   float *deriv, 
	   long nx, 
	   long ny,							\
	   float *xout, 
	   float *yout, 
	   long npt, 
	   long *nvalidx, 
	   float *res)
{
  long i;
  for (i=0;i<=npt;i++) _splin2(xin,yin,im,deriv,nx,ny,nvalidx, \
			       xout[i],yout[i],&res[i]);
}


void 
_spline2grid( float *xin, 
	      float *yin, 
	      float *im, 
	      float *deriv, 
	      long nx, 
	      long ny, 
	      float *xout, 
	      float *yout, 
	      long nxout, 
   	      long nyout, 
	      long *nvalidx,
	      float *res)
/* checked indices. runs good and fast. */
{
  long j,ii,jj;
  float *y2tmp,*yytmp;
  long n;
  long m;

  y2tmp = p_malloc(sizeof(float)*ny);
  yytmp = p_malloc(sizeof(float)*ny);

  for (ii=0;ii<=nxout-1;ii++) {/* loop on out x */

    n=0;

    /* fill Y vector for xout(ii) */
    for (j=0;j<=ny-1;j++) {
      m = nvalidx[j];
      _splint(&xin[n],&im[n],&deriv[n],m,xout[ii],&yytmp[j]);
      n += m;
    }
    /* find second derivative */
    _splinf(yin,yytmp,ny,y2tmp);

    /* find and fill interpolated out Y vector for this xout(ii) */
    for (jj=0;jj<=nyout-1;jj++) {
      _splint(yin,yytmp,y2tmp,ny,yout[jj],&res[jj*nxout+ii]);
    }

  }

  p_free(y2tmp);
  p_free(yytmp);

}

