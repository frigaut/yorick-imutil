/*
 * Author: Francois Rigaut
 * 
 * check file for imutil
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
*/   



require,"imutil.i";
require,"util_fr.i";
require,"random_et.i";


func chk_bilinear(void)
{
  // removed the following lines as it introduced an
  // unnecessary dependancy to yorick-z
  //  im = img_read("crane.jpg");
  //  im = float(im);

//  f = openb("crane.dat");
//  restore,f,im;
//  close,f;
  dim = 512; sdim=32;
  im = gaussdev([2,dim,dim]);
  star = exp(-(dist(sdim)/2.)^2.);
  for (i=1;i<=1000;i++) {
    xs = long(dim/2.+gaussdev()*dim/5.);
    ys = long(dim/2.+gaussdev()*dim/5.);
    if (min(_(xs,ys)) < 1) continue;
    if (max(_(xs,ys)+sdim-1) > dim) continue;
    im(xs:xs+sdim-1,ys:ys+sdim-1) += star/(0.01+random());
  }
  
  //  im = float(jpeg_read("polarbear.jpg"))(sum,,)(,::-1);
  //  palette,"gray.gp";
  
  time = 100;
  // tests of speed:
  xi=1; yi=1; dim=100; nreb=2;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  tic; tmp=bilinear(sim,nreb); eltime=tac()*1000.;
  write,format="bilinear(im,nreb) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",dim,dim,
    nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;
  dim = 512;
  tic; tmp=bilinear(im,nreb); eltime=tac()*1000.;
  write,format="bilinear(im,nreb) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",dim,dim,
    nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;


  xi=1; yi=1; dim=100;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  xy = indices(100);
  tic; tmp=bilinear(sim,xy(,,1),xy(,,2)); eltime=tac()*1000.;
  write,format="bilinear(im,xarray,yarray) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",
    dim,dim,nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;

  //test of equality:
  if (allof(sim == bilinear(sim,1))) {
    write,"Checking that bilinear(im,1) = im ... Yes";
  } else { error,"bilinear(im,1) != im"; }
  
  //displays and check of functionalities:
  xi=50; yi=200; dim=200;
  //  xi=220; yi=300; dim=200;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  //  window,wait=1;
  palette,"earth.gp";

  pli,sim;
  pltitle,"Original";
  write,"Original";
  //pause,time;

  fma; pli,bilinear(sim,2);
  pltitle,"bilinear(im,2)";
  write,"bilinear(im,2)";
  //pause,time;

  fma; pli,bilinear(sim,355,401);
  pltitle,"bilinear(im,355,401)";
  write,"bilinear(im,355,401)";
  //pause,time;

  fma; pli,bilinear(sim,401,355);
  pltitle,"bilinear(im,401,355)";
  write,"bilinear(im,401,355)";
  //pause,time;

  x = indgen(2*dim)/2.;
  y = indgen(2*dim)/2.;
  fma; pli,bilinear(sim,x,y,grid=1);
  pltitle,"bilinear(im,x,y,grid=1)";
  write,"bilinear(im,x,y,grid=1)";
  //pause,time;

  x = indgen(2*dim)/2.;
  y = indgen(2*dim)/3.;
  fma; plg,bilinear(im,x,y);
  pltitle,"bilinear(im,vector!_x,vector!_y)";
  write,"bilinear(im,vector_x,vector_y)";
  //pause,time;

  //speed test:
  //  window,style="nobox.gs",wait=1;
  animate,1;
  for (i=10;i<=dim;i+=4) {fma;pli,bilinear(sim,i,i),cmin=0,cmax=255;}
  for (i=dim/2-1;i>=2;i-=4) {
    fma;
    pli,bilinear(sim(dim/2-i:dim/2+i,dim/2-i:dim/2+i),dim,dim),cmin=0,cmax=255;
  }
  animate,0;
  //pause,time;

  // check of other possibilities for errors:
  //  window,style="work.gs";
  tmp = bilinear(sim,10,300);
  tmp = bilinear(sim,300,10);
  tmp = bilinear(sim,300,300,minus_one=1);
  tmp = bilinear(sim,3,minus_one=1);

  xy = indices(512)-256;
  animate,1;
  for (a=0;a<=180;a+=20) {
    x =  cos(a*pi/180)*xy(,,1) + sin(a*pi/180)*xy(,,2) + 256;
    y = -sin(a*pi/180)*xy(,,1) + cos(a*pi/180)*xy(,,2) + 256;
    fma; pli,bilinear(im,x,y,outside=10);
  }
  animate,0;
  //pause,time;

  
}

func chk_spline2(void)
{
  // removed the following lines as it introduced an
  // unnecessary dependancy to yorick-z
  //  im = img_read("crane.jpg");
  //  im = float(im);

  //f = openb("crane.dat");
  //restore,f,im;
  //close,f;
  dim = 512; sdim=32;
  im = gaussdev([2,dim,dim]);
  star = exp(-(dist(sdim)/2.)^2.);
  for (i=1;i<=1000;i++) {
    xs = long(dim/2.+gaussdev()*dim/5.);
    ys = long(dim/2.+gaussdev()*dim/5.);
    if (min(_(xs,ys)) < 1) continue;
    if (max(_(xs,ys)+sdim-1) > dim) continue;
    im(xs:xs+sdim-1,ys:ys+sdim-1) += star/(0.01+random());
  }

  //  im = float(jpeg_read("polarbear.jpg"))(sum,,)(,::-1);
  //  palette,"gray.gp";
  
  time = 100;
  // tests of speed:
  xi=1; yi=1; dim=100; nreb=2;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  tic; tmp=spline2(sim,nreb); eltime=tac()*1000.;
  write,format="spline2(im,nreb) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",dim,dim,
    nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;
  dim = 512;
  tic; tmp=spline2(im,nreb); eltime=tac()*1000.;
  write,format="spline2(im,nreb) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",dim,dim,
    nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;


  xi=1; yi=1; dim=100;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  xy = indices(100);
  tic; tmp=spline2(sim,xy(,,1),xy(,,2)); eltime=tac()*1000.;
  write,format="spline2(im,xarray,yarray) (%d,%d) -> (%d,%d) : %.1fms, %.2fus/pixel\n",
    dim,dim,nreb*dim,nreb*dim,eltime,eltime*1000/(nreb*dim)^2.;

  //test of equality:
  if (allof(sim == spline2(sim,1))) {
    write,"Checking that spline2(im,1) = im ... Yes";
  } else { error,"spline2(im,1) != im"; }
  
  //displays and check of functionalities:
  xi=50; yi=200; dim=200;
  //  xi=220; yi=300; dim=200;
  sim = im(xi:xi+dim-1,yi:yi+dim-1);
  //  window,wait=1;

  fma;
  pli,sim;
  pltitle,"Original";
  write,"Original";
  //pause,time;

  fma; pli,spline2(sim,2);
  pltitle,"spline2(im,2)";
  write,"spline2(im,2)";
  //pause,time;

  fma; pli,spline2(sim,355,401);
  pltitle,"spline2(im,355,401)";
  write,"spline2(im,355,401)";
  //pause,time;

  fma; pli,spline2(sim,401,355);
  pltitle,"spline2(im,401,355)";
  write,"spline2(im,401,355)";
  //pause,time;

  x = indgen(2*dim)/2.;
  y = indgen(2*dim)/2.;
  fma; pli,spline2(sim,x,y,grid=1);
  pltitle,"spline2(im,x,y,grid=1)";
  write,"spline2(im,x,y,grid=1)";
  //pause,time;
  return;

  x = indgen(2*dim)/2.;
  y = indgen(2*dim)/3.;
  fma; plg,spline2(im,x,y);
  pltitle,"spline2(im,vector!_x,vector!_y)";
  write,"spline2(im,vector_x,vector_y)";
  //pause,time;

  //speed test:
  //  window,style="nobox.gs";
  animate,1;
  for (i=10;i<=dim;i+=4) {fma;pli,spline2(sim,i,i),cmin=0,cmax=255;}
  for (i=dim/2-1;i>=2;i-=4) {
    fma;
    pli,spline2(sim(dim/2-i:dim/2+i,dim/2-i:dim/2+i),dim,dim),cmin=0,cmax=255;
  }
  animate,0;
  //pause,time;

  // check of other possibilities for errors:
  //  window,style="work.gs";
  tmp = spline2(sim,10,300);
  tmp = spline2(sim,300,10);
  tmp = spline2(sim,300,300,minus_one=1);
  tmp = spline2(sim,3,minus_one=1);

  xy = indices(50)-25;
  animate,1;
  for (a=0;a<=180;a+=20) {
    x =  cos(a*pi/180)*xy(,,1) + sin(a*pi/180)*xy(,,2) + 100;
    y = -sin(a*pi/180)*xy(,,1) + cos(a*pi/180)*xy(,,2) + 100;
    fma; pli,spline2(sim,x,y);
  }
  animate,0;
  //pause,time;

  
}

window,wait=1;

tic; g = dist(512); t1=tac();
tic; g = __dist(512); t2=tac();
write,format="Create dist(512): %fs (with interpreted function %f)\n",t1,t2;

tic; gc = clip(g,10,50); t1=tac();
tic; gc = min(max(g,10),50); t2=tac();
write,format="clip image 512x512: %fs (with Buildin yorick function %f)\n",t1,t2;

tic; gc = eclat(g); t1=tac();
tic; gc = roll(g); t2=tac();
write,format="Roll image 512x512: %fs (with Buildin yorick function %f)\n",t1,t2;

tic; greb = bin2(g); t=tac();
write,format="rebin 512x512 -> 256x256: %fs\n",t;

tic; x = gaussdev(100000); t1=tac();
tic; x = random_n(100000); t2=tac();
write,format="100000 normal random numbers: %fs (interpreted function %fs)\n",t1,t2;

v = array(10,100000);
tic; x = poidev(v); t1=tac();
tic; x = random_poisson(v); t2=tac();
write,format="Poisson of array(10,100000): %fs (interpreted function %fs)\n",t1,t2;

//Interpolation tests

chk_bilinear;
chk_spline2;
write,"Image interpolation tests OK";

// Sedgesort tests:

n = 500000;
write,format="Sorting %d numbers with sedgesort\n",n;

a = long(10000*random(n));
tic; v = a(sort(a));   t1=tac();
tic; v = sedgesort(a); t2=tac();
write,format="Long:   Regular sort: %fs / sedgesort: %fs\n",t1,t2;

a = float(10000*random(n));
tic; v = a(sort(a));   t1=tac();
tic; v = sedgesort(a); t2=tac();
write,format="Float:  Regular sort: %fs / sedgesort: %fs\n",t1,t2;

a = double(10000*random(n));
tic; v = a(sort(a));   t1=tac();
tic; v = sedgesort(a); t2=tac();
write,format="Double: Regular sort: %fs / sedgesort: %fs\n",t1,t2;

write,"All tests OK";

