/* imutil Yorick plugin
 * $Id$
 * Francois Rigaut, 2003-2011
 *
 * Copyright (c) 2003, Francois RIGAUT (frigaut@gemini.edu, Gemini
 * Observatory, 670 N A'Ohoku Place, HILO HI-96720).
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
 * $log$
 *
*/

plug_in,"imutil";
require,"util_fr.i";

//func is_scalar(x) { return (is_array(x) && ! dimsof(x)(1)); }
//func is_vector(x) { return (is_array(x) && dimsof(x)(1) == 1); }

func dpc(image,badpixmap,silent=)
/* DOCUMENT dpc(image,badpixmap,silent=)
  *Correction of bad pixels in an image by averaging the (good) neighbors.
 * Dead pixel correction.
 * image     : image to process
 * badpixmap : bad pixel map (same dimension as image). 1 means bad pixel.
 * silent      be silent. Default is to write out # of pixels to process.
 * F.Rigaut 2011/03/16
 * SEE ALSO: deadpix() from yorick-yutils. This function superseeds deadpix()
 * and is several order of magnitudes faster.
 */
{
  // the bad pixels we replace are going to be approximate, so float is enough
  ima    = float(image);
  gpm    = short(1-badpixmap);
  w0     = where(badpixmap);
  if (numberof(w0)==0) return image;
  w      = w0-1; // -1 = yorick 1-based vs C 0-based
  w      = long(w);
  w2     = w*0;
  stride = long(dimsof(ima)(2));
  npix   = long(numberof(ima));
  nw     = long(numberof(w));
  if (!silent) write,format="dpc: %d dead pixels to process...",nw;
  tic,9;
  npassmax = 100l;
  _dpc,&ima,&gpm,&w,&w2,stride,npix,nw,npassmax;
  image(w0) = ima(w0); // will be cast automatically
  if (!silent) write,format="Done (%.3fs)\n",tac(9);
  return image;
}
if (is_func(deadpix)) deadpix=dpc;


func sigfil(image,nsigma,iter=,silent=)
/* DOCUMENT func sigfil(image,nsigma,iter=,silent=)
   Filter out the pixels that deviate from the local statistics.
   The mean and rms of the 8 (the minimum and maximum of these
   8 neighbors are excluded in the mean and rms computation) is
   computed. All pixels that deviates more than "nsigma" rms
   from the mean are flagged as bad pixels. The image and newly
   created bad pixel map are passed to the routine dpc()
   for correction. The processus can be iterated.
   image	: input image
   nsigma: number of rms about the local mean out of which is
   pixel is considered aberrant. nsigma >= 5 recommended.
   iter  : Keyword, number of iterations. Recommended value : 3-5
   silent: No verbose
   F.Rigaut 2011/03/16
   This function superseeds sigmaFilter from yorick-yutils
   SEE ALSO: dpc
  */
{
  if (!iter) iter=100;
  ima    = float(image);
  stride = long(dimsof(ima)(2));
  npix   = long(numberof(ima));
  abpm   = array(short,dimsof(ima));
  do {
    bpm    = array(short,dimsof(ima));
    _sigfil,&ima,&bpm,float(nsigma),stride,npix;
    abpm += bpm;
    nbp = sum(bpm);
    ima = dpc(ima,bpm,silent=silent);
    iter--;
  } while ((nbp)&&(iter>0));
  w = where(abpm);
  if (numberof(w)) image(w) = ima(w);
  return image;
}
if (is_func(sigmaFilter)) sigmaFilter=sigfil;


func cart2pol(image,&r,&theta,xc=,yc=,ntheta=,nr=,tor=,splin=,outside=)
/* DOCUMENT cart2pol(image,&r,&theta,xc=,yc=,ntheta=,nr=,tor=,
                     splin=,outside=)

   Cartesian to Polar coordinate coordinate mapping

   image : input image, in cartesian coorindate (square grid)

   optional output: r: 1D vector to hold the r coordinates at
   which output image is mapped theta: same for theta

   KEYWORDS:
   xc, yc: Center for coordinate transform. Note that this is
   compatible with the center defined by dist(), or rotate2(), but is
   offset by 0.5 pixels w.r.t what you read on the yorick graphic
   window. I.e. the center of the bottom- left pixel is (1,1) in this
   function's conventions, not (0.5,0.5).

   nr, ntheta: number of points over the output radius and theta range

   tor: upper output radius value (in pixels)

   splin: use spline2() instead of bilinear() for the interpolation

   outside: value for outliers.

   SEE ALSO: rotate2, spline2, bilinear
*/
{
  d = dimsof(image);
  if (xc==[]) xc=ceil(d(2)/2.+0.5);
  if (yc==[]) yc=ceil(d(3)/2.+0.5);
  if (!tor) {
    tor=max(abs(d(2:3)-_(xc,yc)));
  }
  if (!ntheta) ntheta=max(d(2:3));
  if (!nr) nr = long(ceil(tor)+1); else nr=long(nr);
  r     = array(1.,ntheta)(-,)*span(0.,tor,nr);
  theta = span(0.,2*pi,ntheta+1)(1:-1)(-,)*array(1.,nr);

  x = r*cos(theta)+xc;
  y = r*sin(theta)+yc;

  r = r(,1);
  theta = theta(1,);

  if (splin) return spline2(image,x,y,outside=outside);
  return bilinear(image,x,y,outside=outside);
}



func rotate2(image,angle,xc=,yc=,splin=,outside=)
/* DOCUMENT rotate2(image,angle,xc=,yc=,splin=,outside=)

   Rotate the input image. Angle is in degrees, CCW.

   KEYWORDS:
   xc, yc: Center for coordinate transform. Note that this is
   compatible with the center defined by dist(), or cart2pol(), but is
   offset by 0.5 pixels w.r.t what you read on the yorick graphic
   window. I.e. the center of the bottom- left pixel is (1,1) in this
   function's conventions, not (0.5,0.5).

   splin: use spline2() instead of bilinear() for the interpolation

   outside: value for outliers.

   SEE ALSO: spline2, cart2pol, bilinear
*/
{
  angle *= pi/180.;

  d = dimsof(image);

  xy = indices(d);

  if (xc==[]) xc=ceil(d(2)/2.+0.5);
  if (yc==[]) yc=ceil(d(3)/2.+0.5);

  xy(,,1)-=xc;
  xy(,,2)-=yc;

  x =  cos(angle)*xy(,,1) + sin(angle)*xy(,,2);
  y = -sin(angle)*xy(,,1) + cos(angle)*xy(,,2);

  x +=xc;
  y +=yc;

  if (splin) return spline2(image,x,y,outside=outside);
  return bilinear(image,x,y,outside=outside);
}



func bin2d(in,binfact)
/* DOCUMENT func bin2d(in,binfact)
   Returns the input 2D array "in", binned with the binning factor "binfact".
   The input array X and/or Y dimensions needs not to be a multiple of
   "binfact"; The final/edge pixels are in effect replicated if needed.
   This routine prepares the parameters and calls the C routine _bin2d.
   The input array can be of type long, float or double.
   Last modified: Dec 15, 2003.
   Author: F.Rigaut
   SEE ALSO: _bin2d
 */

{
  if ( (binfact-long(binfact)) != 0) write,"*** Warning: binfact has to be an int";
  if (binfact<1) error,"binfact has to be >= 1";

  binfact = int(binfact); // this *has* to be a int.

  nx = int(dimsof(in)(2));
  ny = int(dimsof(in)(3));
  fx = int(ceil(nx/float(binfact)));
  fy = int(ceil(ny/float(binfact)));

  // Test type of array and call appropriate routine
  if (typeof(in) == "int") in=long(in);

  if (typeof(in) == "long") {

    // define/allocate output image
    outData = array(long,[2,fx,fy]);
    err = _bin2d_long(&in,nx,ny,&outData,fx,fy,binfact);
    return outData;

  } else if (typeof(in) == "float") {

    // define/allocate output image
    outData = array(float,[2,fx,fy]);
    err = _bin2d_float(&in,nx,ny,&outData,fx,fy,binfact);
    return outData;

  } else if (typeof(in) == "double") {

    // define/allocate output image
    outData = array(double,[2,fx,fy]);
    err = _bin2d_double(&in,nx,ny,&outData,fx,fy,binfact);
    return outData;

  } else { error,"Unsupported Data Type"; }

}




func dist(dim,xc=,yc=)
/* DOCUMENT func dist(dim,xc=,yc=)
 * Return an array which elements are the distance to (xc,yc). xc and
 * yc can be omitted, in which case they are defaulted to size/2+1.
 * F.Rigaut, 2003/12/10.
 * SEE ALSO:
*/

{
  dim = long(dim);

  if (is_scalar(dim)) dim=[2,dim,dim];
  if ((is_vector(dim)) && (dim(1)!=2))
    error,"Dist only deals with 2D square arrays";

  d = array(float,dim);

  if (xc!=[]) {xc = float(xc-1.);} else {xc = float(dim(2)/2);}
  if (yc!=[]) {yc = float(yc-1.);} else {yc = float(dim(3)/2);}

  res = _dist(&d,dim(2),dim(3),xc,yc);

  return d;
}


func clip(inarray,xmin,xmax)
/* DOCUMENT func clip(inarray, mini, maxi);
 * Returns the argument, which has been "clipped" to mini
 * and maxi, i.e. in which all elements lower than "mini"
 * have been replaced by "mini" and all elements greater
 * than "maxi" by "maxi".
 * Either "mini" and "maxi" can be ommited, in which case
 * the corresponding mini or maxi is not clipped.
 * Equivalent to the IDL ">" and "<" operators.
 *
 * Can clip in place
 * prompt> clip,a,0.,1.
 * or out of place
 * prompt> res = clip(a,0.,1.)
 *
 * F.Rigaut, 2001/11/10.
 * 2003/12/4: Using now a C version which is called from this routine.
*/

{
  local x;

  // Check data type to clip. Only supports long, float and double.
  if ( (typeof(inarray) != "char")  & (typeof(inarray) != "short") &
       (typeof(inarray) != "int")   & (typeof(inarray) != "long") &
       (typeof(inarray) != "float") & (typeof(inarray) != "double") ) {
    error,"Unknown Data type in clip";
  }

  sub = am_subroutine();
  if (sub) {eq_nocopy,x,inarray;} else {x = inarray;}

  // Simple check that input limits make sense
  if ( (xmin != []) && (xmax != []) && (xmin > xmax) ) {
    error,"xmin > xmax in clip";
  }

  // Have to deal with missing min or max in this lenghthy way:
  // no min specified:
  if (xmin == []) {
    if (typeof(x) == "char") {
      res = clipmaxchar(&x,char(xmax),numberof(x));
    }
    if (typeof(x) == "short") {
      res = clipmaxshort(&x,short(xmax),numberof(x));
    }
    if (typeof(x) == "int") {
      res = clipmaxint(&x,int(xmax),numberof(x));
    }
    if (typeof(x) == "long") {
      res = clipmaxlong(&x,long(xmax),numberof(x));
    }
    if (typeof(x) == "float") {
      res = clipmaxfloat(&x,float(xmax),numberof(x));
    }
    if (typeof(x) == "double") {
      res = clipmaxdouble(&x,double(xmax),numberof(x));
    }
    return x;
  }

  // no max specified:
  if (xmax == []) {
    if (typeof(x) == "char") {
      res = clipminchar(&x,char(xmin),numberof(x));
    }
    if (typeof(x) == "short") {
      res = clipminshort(&x,short(xmin),numberof(x));
    }
    if (typeof(x) == "int") {
      res = clipminint(&x,int(xmin),numberof(x));
    }
    if (typeof(x) == "long") {
      res = clipminlong(&x,long(xmin),numberof(x));
    }
    if (typeof(x) == "float") {
      res = clipminfloat(&x,float(xmin),numberof(x));
    }
    if (typeof(x) == "double") {
      res = clipmindouble(&x,double(xmin),numberof(x));
    }
    return x;
  }

  // min and max specified:
  if (typeof(x) == "char") {
    res = clipchar(&x,char(xmin),char(xmax),numberof(x));
  }
  if (typeof(x) == "short") {
    res = clipshort(&x,short(xmin),short(xmax),numberof(x));
  }
  if (typeof(x) == "int") {
    res = clipint(&x,int(xmin),int(xmax),numberof(x));
  }
  if (typeof(x) == "long") {
    res = cliplong(&x,long(xmin),long(xmax),numberof(x));
  }
  if (typeof(x) == "float") {
    res = clipfloat(&x,float(xmin),float(xmax),numberof(x));
  }
  if (typeof(x) == "double") {
    res = clipdouble(&x,double(xmin),double(xmax),numberof(x));
  }

  return x;
}



func eclat(image)
/* DOCUMENT func eclat(image)
 * Equivalent, but significantly faster than roll. Transpose the four main
 * quadrants of a 2D array. Mostly used for FFT applications.
 * The C function can be called directly in time critical loops as
 * _eclat_type,&image,nx,ny
 * with type = long, float or double (e.g. _eclat_float,...)
 *
 *  Can invoque in place
 *  prompt> eclat,image
 *  or out of place
 *  prompt> res = eclat(image)
 *
 * F.Rigaut, 2001/11/10.
 * SEE ALSO: roll.
 */

{
  local x;

  nx = int(dimsof(image)(2));
  ny = int(dimsof(image)(3));  //fixed 2->3 may18, 2005.


  sub = am_subroutine();
  if (sub) {eq_nocopy,x,image;} else {x = image;}

  if (typeof(image) == "long") {
    _eclat_long,&x,nx,ny;
  } else if (typeof(image) == "float") {
    _eclat_float,&x,nx,ny;
  } else if (typeof(image) == "double") {
    _eclat_double,&x,nx,ny;
  } else { error,"Unsupported data type";}

  return x;
}



func poidev(vec)
/* DOCUMENT func poidev(vec)
   Compute random values following a Poisson Distribution.
   Input: array containing the desired mean value(s)
   output: randomized array of same dimension
   EXAMPLE:
   p = poidev(array(20.,[2,128,128]));
   returns a 128x128 array with pixels following a Poisson distrib.
   of mean=20.
   SEE ALSO: gaussdev
 */
{
  vec = float(vec);
  _poidev,vec,numberof(vec);
  return vec;
}



func gaussdev(dims)
/* DOCUMENT func gaussdev(dims)
   Returns an array (as specified by dims) of random values
   following a normal distribution of mean=0 and standard
   deviation = 1.

   EXAMPLES:

   gaussdev(100) returns a 100 element vector
   gaussdev([2,512,512]) returns a 512x512 array

   print,gaussdev([2,512,512])(rms)
   0.991666

   SEE ALSO: poidev
 */
{
  vec = array(float,dims);
  _gaussdev,vec,numberof(vec);
  return vec;
}



func spline2(image,arg1,arg2,grid=,minus_one=,outside=,mask=)
/* DOCUMENT spline2(image,arg1,arg2,grid=,minus_one=,outside=,mask=)

  Interpolate regularly sampled 2D arrays using 2D spline.

        spline2(image,nrebin)
    or  spline2(image,dimx,dimy)
    or  spline2(image,x,y,grid=1)
    or  spline2(image,x,y)

  spline2(image,nrebin)
    Returns image interpolated by a factor N
    Output dimension = N * input dimension
    The whole image is interpolated

  spline2(image,dimx,dimy)
    Returns interpolated image
    Output dimension = dimx * dimy
    The whole image is interpolated

  spline2(image,x,y,grix=1)
    Returns interpolated image on a grid of points specified
    by the coordinates x & y. The output is a 2D array
    of dimension numberof(x) * numberof(y).

  spline2(image,x,y)
    Returns image interpolated at points of coordinates (x,y).
    X and Y can be of any dimension, e.g. scalar (but in that
    case it should be entered as a one element vector to avoid
    confusion with the second form, e.g. spline2(im,[245.34],[45.3])),
    vectors (for instance along a line of circle) or can be 2D
    array to specify an arbitrary transform (e.g. a rotation).
    This is the general form of spline2. It is significantly
    slower than the previous forms that take advantage of the
    cartesian nature of the output coordinates, so don't
    use it unless it is necessary.

  NOTE on the input coordinates: In the input image, the lower
  left pixel has coordinates [xin(1),yin(1)]=[1,1]. Pixels are
  spaced by one unit, so [xin(2),yin(1)]=[2,1],...

  KEYWORD minus_one specify that the last column/row should not
  be extrapolated, i.e. the output image extend from xin(1) to
  xin(0) inclusive, not beyond (same for y).

  KEYWORD mask can be used if some of the input image pixels
  are invalid and not to be used to compute the output image.
  Bad pixel can be interpolated using mask. Mask is an array
  of same dimension as image. 0 mark an invalid data point.

  EXAMPLES: if "image" is a 512x512 image:

  spline2(image,2) returns image interpolated on a 1024x1024 grid

  spline2(image,1024,800) returns image interpolated on a 1024x800 grid

  xreb=100+indgen(400)/2.;
  yreb=50.4+indgen(200)/2.;
  spline2(image,xreb,yreb,grid=1) returns image, interpolated in 2D
  over a XY grid defined by xreb and yreb.

  spline2(image,[300.3],[200.4]) returns the interpolation of image
  at point [300.3,200.4]

  spline2(image,indgen(512),indgen(512)) returns a vector of values
  interpolated at [x,y] = [[1,1],[2,2],...,[512,512]]

  x = (indgen(400)/3.+50)*array(1,200)(-,);
  y = (array(1,400)*(indgen(200)/2.+55)(-,));
  spline2(image,x,y) returns a 2D interpolated array at indices defined
  by x and Y

  xy = indices(512);
  alpha = 34*pi/180.;
  x = cos(alpha)*xy(,,1)+sin(alpha)*xy(,,2);
  y = -sin(alpha)*xy(,,1)+cos(alpha)*xy(,,2);
  spline2(image,x,y) returns image rotated by alpha
  SEE ALSO: bilinear, spline
 */
{
  if (structof(image) != float) image=float(image);

  d = dimsof(image);
  nx = d(2);
  ny = d(3);

  xy  = indices([2,nx,ny]);
  xin = float(xy(,,1));
  yin = float(indgen(ny));
  xy = [];

  if (is_void(mask)) mask = short(image*0+1);
  else mask = short(mask);

  deriv = image*0.f;

  if (is_void(arg2)) {

    //case 1: spline2(image,nrebin)
    nreb = long(arg1);
    xreb = float((indgen(nx*nreb)-1.f)/nreb+1.0);
    yreb = float((indgen(ny*nreb)-1.f)/nreb+1.0);

  } else if (is_scalar(arg1) && is_scalar(arg2)) {

    //case 2: spline2(image,dimx,dimy)
    dimx = long(arg1);
    dimy = long(arg2);
    xreb = float((indgen(dimx)-1.f)/dimx*nx+1.0);
    yreb = float((indgen(dimy)-1.f)/dimy*ny+1.0);

  } else if (grid == 1) {

    // case 3: spline2(image,x,y,grid=1)
    xreb = float(arg1);
    yreb = float(arg2);

  } else {

    // gotta be case 4: spline2(image,x,y)
    xreb = float(arg1);
    yreb = float(arg2);
    if (anyof(dimsof(xreb) != dimsof(yreb))) {
      error,"X and Y have to have the same dimension";
    }
    npt = numberof(xreb);
    _s2gen = 1;

  }

  if (minus_one) {
    xreb = xreb(where(xreb <= xin(0)));
    yreb = yreb(where(yreb <= yin(0)));
  }
  nxreb = numberof(xreb);
  nyreb = numberof(yreb);

  nvalidx = long(mask(sum,));
  wm = where(mask(*));
  xin = xin(wm);
  image = image(wm);
  wm = [];

  _splie2,xin,image,nx,ny,deriv,nvalidx;

  if (_s2gen) {

    res = array(float,dimsof(xreb));
    _spline2,xin,yin,image,deriv,nx,ny,xreb,yreb,npt,nvalidx,res;

    if (!is_void(outside)) {
      _mask = where((xreb>xin(0))|(xreb<xin(1))|(yreb>yin(0))|(yreb<yin(1)));
      if (numberof(_mask)) res(_mask) = outside;
    }

  } else {

    res = array(float,[2,nxreb,nyreb]);
    _spline2grid,xin,yin,image,deriv,nx,ny,xreb,yreb,nxreb,nyreb,nvalidx,res;

    if (!is_void(outside)) {
      maskx = where((xreb>xin(0))|(xreb<xin(1)));
      if (numberof(maskx)) res(maskx,) = outside;
      masky = where((yreb>yin(0))|(yreb<yin(1)));
      if (numberof(masky)) res(,masky) = outside;
    }

  }

  return res;
}


func bilinear(image,arg1,arg2,grid=,minus_one=,outside=)
/* DOCUMENT bilinear(image,arg1,arg2,grid=,minus_one=,outside=)

  Interpolate regularly sampled 2D arrays using bilinear interpolation.

        bilinear(image,nrebin)
    or  bilinear(image,dimx,dimy)
    or  bilinear(image,x,y,grid=1)
    or  bilinear(image,x,y)

  bilinear(image,nrebin)
    Returns image interpolated by a factor N
    Output dimension = N * input dimension
    The whole image is interpolated

  bilinear(image,dimx,dimy)
    Returns interpolated image
    Output dimension = dimx * dimy
    The whole image is interpolated

  bilinear(image,x,y,grix=1)
    Returns interpolated image on a grid of points specified
    by the coordinates x & y. The output is a 2D array
    of dimension numberof(x) * numberof(y).

  bilinear(image,x,y)
    Returns image interpolated at points of coordinates (x,y).
    X and Y can be of any dimension, e.g. scalar (but in that
    case it should be entered as a one element vector to avoid
    confusion with the second form, e.g. bilinear(im,[245.34],[45.3])),
    vectors (for instance along a line of circle) or can be 2D
    array to specify an arbitrary transform (e.g. a rotation).
    This is the general form of bilinear.

  NOTE on the input coordinates: In the input image, the lower
  left pixel has coordinates [xin(1),yin(1)]=[1,1]. Pixels are
  spaced by one unit, so [xin(2),yin(1)]=[2,1],...

  KEYWORD minus_one specify that the last column/row should not
  be extrapolated, i.e. the output image extend from xin(1) to
  xin(0) inclusive, not beyond (same for y).

  EXAMPLES: if "image" is a 512x512 image:

  bilinear(image,2) returns image interpolated on a 1024x1024 grid

  bilinear(image,1024,800) returns image interpolated on a 1024x800 grid

  xreb=100+indgen(400)/2.;
  yreb=50.4+indgen(200)/2.;
  bilinear(image,xreb,yreb,grid=1) returns image, interpolated in 2D
  over a XY grid defined by xreb and yreb.

  bilinear(image,[300.3],[200.4]) returns the interpolation of image
  at point [300.3,200.4]

  bilinear(image,indgen(512),indgen(512)) returns a vector of values
  interpolated at [x,y] = [[1,1],[2,2],...,[512,512]]

  x = (indgen(400)/3.+50)*array(1,200)(-,);
  y = (array(1,400)*(indgen(200)/2.+55)(-,));
  bilinear(image,x,y) returns a 2D interpolated array at indices defined
  by x and Y

  xy = indices(512);
  alpha = 34*pi/180.;
  x = cos(alpha)*xy(,,1)+sin(alpha)*xy(,,2);
  y = -sin(alpha)*xy(,,1)+cos(alpha)*xy(,,2);
  bilinear(image,x,y) returns image rotated by alpha
  SEE ALSO: spline2, interp2
*/
{
  if (structof(image) != float) image=float(image);

  if (is_void(outside)) {
    skipoutside=0l;
    outside=0.f;
  } else {
    skipoutside=1l;
    outside = float(outside);
  }

  d = dimsof(image);
  if (d(1)!=2) error,"2D array only";
  nx = d(2);
  ny = d(3);
  xin = indgen(nx);
  yin = indgen(ny);

  if (is_void(arg2)) {  //case 1: nreb
    nreb = long(arg1);
    xreb = float((indgen(nx*nreb)-1.f)/nreb+1.0);
    yreb = float((indgen(ny*nreb)-1.f)/nreb+1.0);
  } else if (is_scalar(arg1) && is_scalar(arg2)) { //case 2: dimx,dimy
    dimx = long(arg1);
    dimy = long(arg2);
    xreb = float((indgen(dimx)-1.f)/dimx*nx+1.0);
    yreb = float((indgen(dimy)-1.f)/dimy*ny+1.0);  // bug detected and corrected 2007jun14 (nx->ny)
  } else if (grid == 1) { // case 3
    xreb = float(arg1);
    yreb = float(arg2);
  } else { // gotta be case 4
    xreb = float(arg1);
    yreb = float(arg2);
    if (anyof(dimsof(xreb) != dimsof(yreb))) {
      error,"X and Y have to have the same dimension";
    }
    npt = numberof(xreb);
    _i2gen = 1;
  }

  if (minus_one) {
    xreb = xreb(where(xreb <= xin(0)));
    yreb = yreb(where(yreb <= yin(0)));
  }

  if (!_i2gen) {
    nxreb = numberof(xreb);
    nyreb = numberof(yreb);
    xreb = xreb*array(1.f,nyreb)(-,);
    yreb = array(1.f,nxreb)*yreb(-,);
  }

  res = array(outside,dimsof(xreb));

  _bilinear,image,nx,ny,res,xreb,yreb,numberof(xreb),skipoutside;

  return res;
}


func sedgesort(vec)
/* DOCUMENT func sedgesort(vec)
   Fast sort method. Typically 2 to 5 times faster than the yorick
   stock function. Only valid input is a 1D array.
   WARNING! Returns the sorted vector. This is different than the
   stock yorick sort, which returns a vector of indices.
   SEE ALSO: sort, sedgemedian
 */
{
  local x;
  if (dimsof(vec)(1) != 1) {
    error,"sedgesort only works on 1d array";
  }
  if (noneof(typeof(vec)==["short","long","float","double"])) {
    error,swrite(format="sorry, type %s not supported",typeof(x));
  }
  x = _(vec,max(vec)+1);
  xmin = min(x);
  x = x - xmin;
  if      (typeof(x)=="short")  {_sedgesort_short ,&x,int(numberof(x));}
  else if (typeof(x)=="long")   {_sedgesort_long  ,&x,int(numberof(x));}
  else if (typeof(x)=="float")  {_sedgesort_float ,&x,int(numberof(x));}
  else if (typeof(x)=="double") {_sedgesort_double,&x,int(numberof(x));}
  return x(1:-1)+xmin;
}


func sedgemedian(vec)
/* DOCUMENT func sedgemedian(vec)
   Returns the median of the input array (1D). Uses sedgesort fast
   sort.
   SEE ALSO: sedgesort, median
 */
{
  if (dimsof(vec)(1) != 1) {
    error,"sedgesort only works on 1d array";
  }
  if (noneof(typeof(vec)==["short","long","float","double"])) {
    error,swrite(format="sorry, type %s not supported",typeof(x));
  }
  x = _(vec,max(vec)+1);
  xmin = min(x);
  x = x - xmin;
  if      (typeof(x)=="short")  {_sedgesort_short ,&x,int(numberof(x));}
  else if (typeof(x)=="long")   {_sedgesort_long  ,&x,int(numberof(x));}
  else if (typeof(x)=="float")  {_sedgesort_float ,&x,int(numberof(x));}
  else if (typeof(x)=="double") {_sedgesort_double,&x,int(numberof(x));}
  x = x(1:-1);
  len = numberof(x);
  if (len%2) {
    //odd# of elements.
    return x(len/2+1)+xmin;
  } else {
    //even# of elements. must return avg value of middle elements
    return (x(len/2)+x(len/2+1))/2.f+xmin;
  }
}



func cpc(im,fmin,fmax)
/* DOCUMENT func cpc(im,fmin,fmax)
   return clipped image at
   from cmin=fraction fmin of pixel intensity
   to   cmax=fraction fmax of pixel intensity
   0 <= fmin < fmax <= 1
   example: pli,cpc(im)
   SEE ALSO:
 */
{
  s = sedgesort(im(*));
  n = numberof(im);
  if (fmin==[]) fmin = 0.10;
  if (fmax==[]) fmax = 0.995;
  if ((fmin<0)||(fmax<0)) error,"fmin and fmax should be > 0"
  if ((fmin>1)||(fmax>1)) error,"fmin and fmax should be < 1"
  if (fmin>fmax) error,"fmin should be < fmax"
  x1=s(long(round(n*fmin)));
  x2=s(long(round(n*fmax)));
  return clip(im,x1,x2);
}


extern _dpc
/* PROTOTYPE
   void _dpc(pointer im, pointer gpm, pointer w, pointer w2, long stride, long npix, long nw, long npassmax)
*/

extern _sigfil
/* PROTOTYPE
   void _sigfil(pointer im, pointer bpm, float nsigma, long stride, long npix)
*/

extern _bin2d_long
/* PROTOTYPE
   int _bin2d_long(pointer in, int nx, int ny, pointer out, int fx, int fy, int binfact)
*/

extern _bin2d_float
/* PROTOTYPE
   int _bin2d_float(pointer in, int nx, int ny, pointer out, int fx, int fy, int binfact)
*/

extern _bin2d_double
/* PROTOTYPE
   int _bin2d_double(pointer in, int nx, int ny, pointer out, int fx, int fy, int binfact)
*/

extern _dist;
/* PROTOTYPE
   void _dist(pointer dptr, long dimx, long dimy, float xc, float yc)
*/

extern clipchar;
/* PROTOTYPE
   int clipchar(pointer x, char xmin, char xmax, long n)
*/

extern clipshort;
/* PROTOTYPE
   int clipshort(pointer x, short xmin, short xmax, long n)
*/

extern clipint;
/* PROTOTYPE
   int clipint(pointer x, int xmin, int xmax, long n)
*/

extern cliplong;
/* PROTOTYPE
   int cliplong(pointer x, long xmin, long xmax, long n)
*/

extern clipfloat;
/* PROTOTYPE
   int clipfloat(pointer x, float xmin, float xmax, long n)
*/

extern clipdouble;
/* PROTOTYPE
   int clipdouble(pointer x, double xmin, double xmax, long n)
*/

extern clipminchar;
/* PROTOTYPE
   int clipminchar(pointer x, char xmin, long n)
*/

extern clipminshort;
/* PROTOTYPE
   int clipminshort(pointer x, short xmin, long n)
*/

extern clipminint;
/* PROTOTYPE
   int clipminint(pointer x, int xmin, long n)
*/

extern clipminlong;
/* PROTOTYPE
   int clipminlong(pointer x, long xmin, long n)
*/

extern clipminfloat;
/* PROTOTYPE
   int clipminfloat(pointer x, float xmin, long n)
*/

extern clipmindouble;
/* PROTOTYPE
   int clipmindouble(pointer x, double xmin, long n)
*/

extern clipmaxchar;
/* PROTOTYPE
   int clipmaxchar(pointer x, char xmax, long n)
*/

extern clipmaxshort;
/* PROTOTYPE
   int clipmaxshort(pointer x, short xmax, long n)
*/

extern clipmaxint;
/* PROTOTYPE
   int clipmaxint(pointer x, int xmax, long n)
*/

extern clipmaxlong;
/* PROTOTYPE
   int clipmaxlong(pointer x, long xmax, long n)
*/

extern clipmaxfloat;
/* PROTOTYPE
   int clipmaxfloat(pointer x, float xmax, long n)
*/

extern clipmaxdouble;
/* PROTOTYPE
   int clipmaxdouble(pointer x, double xmax, long n)
*/

extern _eclat_long
/* PROTOTYPE
   void _eclat_long(pointer ar, int nx, int ny)
*/

extern _eclat_float
/* PROTOTYPE
   void _eclat_float(pointer ar, int nx, int ny)
*/

extern _eclat_double
/* PROTOTYPE
   void _eclat_double(pointer ar, int nx, int ny)
*/

extern _mynoop1
/* PROTOTYPE
   int _mynoop1(void)
*/

extern _gaussdev
/* PROTOTYPE
   void _gaussdev(float array xm, long n)
*/

extern _poidev
/* PROTOTYPE
   void _poidev(float array xm, long n)
*/

extern ran1init
/* PROTOTYPE
   void ran1init(void)
*/

extern _spline2grid
/* PROTOTYPE
  void _spline2grid(float array xin, float array yin, float array image,
  float array deriv, long nx, long ny,
  float array xreb, float array yreb, long nxreb, long nyreb,
  long array nvalidx, float array res)
*/

extern _spline2
/* PROTOTYPE
  void _spline2(float array xin, float array yin, float array image,
  float array deriv, long nx, long ny,
  float array xout, float array yout, long npt, long array nvalidx,
  float array res)
*/

extern _splie2
/* PROTOTYPE
   void _splie2(float array x, float array image,
   long nx, long ny, float array deriv, long array nvalidx)
*/

extern _bilinear
/* PROTOTYPE
   void _bilinear(float array image, long nx, long ny, float array out,
   float array xout, float array yout, long nout, long skipoutside)
*/

extern _sedgesort_long
/* PROTOTYPE
  void _sedgesort_long(pointer dataprt, int len)
*/
extern _sedgesort_float
/* PROTOTYPE
  void _sedgesort_float(pointer dataprt, int len)
*/
extern _sedgesort_double
/* PROTOTYPE
  void _sedgesort_double(pointer dataptr, int len)
*/
extern _sedgesort_short
/* PROTOTYPE
  void _sedgesort_short(pointer dataptr, int len)
*/

// comment following line to have a deterministic random (!) start...
ran1init;  // init random function for poidev.

