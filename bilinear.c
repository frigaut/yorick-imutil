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


void
_bilinear(float *image, long nx, long ny,
	  float *out, float *xout, float *yout, long nout, long skipoutside)

{
  long i;
  long i0,j0,i1,j1,i00,i01,i10,i11;
  float wi,wj,w00,w01,w10,w11;

    /* Loop on indices of output image */
  for (i=0;i<nout;i++) {

    if (((xout[i]<1)|(xout[i]>nx)|(yout[i]<1)|(yout[i]>ny))&skipoutside) continue;

    i0 = (long)(xout[i])-1; /* -1 because C indices are zero-based */
    j0 = (long)(yout[i])-1;
    i1 = i0+1;
    j1 = j0+1;

    if (i0<0)      i0=0;
    if (i0>(nx-1)) i0=nx-1;
    if (j0<0)      j0=0;
    if (j0>(ny-1)) j0=ny-1;

    if (i1<0)      i1=0;
    if (i1>(nx-1)) i1=nx-1;
    if (j1<0)      j1=0;
    if (j1>(ny-1)) j1=ny-1;

    /* global index = col# + Ncolumns * row# */
    i00 = i0+nx*j0;
    i10 = i1+nx*j0;
    i01 = i0+nx*j1;
    i11 = i1+nx*j1;

    /* Computes the weights for the 4 surrounding pixels */
    wi = 1.0f-(xout[i]-(long)(xout[i]));
    wj = 1.0f-(yout[i]-(long)(yout[i]));

    w00 = wi*wj;
    w10 = (1-wi)*wj;
    w01 = wi*(1-wj);
    w11 = (1-wi)*(1-wj);

    /* Finaly, compute and integrate outphase */
    out[i] = image[i00]*w00+image[i10]*w10+image[i01]*w01+image[i11]*w11;
  }
}
