
/*
   FPToolkit.c : A simple set of graphical tools.
   FPToolkitDemo.c 
   Copyright (C) 2018  Ely

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License (version 3)
   as published by the Free Software Foundation.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
   */





/*

   This code assumes an underlying X11 environment.

   Most freshly installed unbuntu environments do not have
   the X11 developr stuff they'll need to do graphics.
   To download X11 developer stuff, connect to the internet and 
   issue the following two commands.  Each will ask for your password
   and each will take a few minutes.  At some point it might even
   look like nothing is happening....be patient :

   sudo  apt-get  install  libx11-dev     

   sudo  apt-get  install  xorg-dev

*/




/*
   If this file and the file, FPToolkit.c, are in the
   same directory that you are trying to compile in,  
   do the following :

   cc  FPToolkitDemo.c   -lm  -lX11

*/


#include  "FPToolkit.c"


int main()
{
  int    swidth, sheight ;
  double lowleftx, lowlefty, width, height ;
  double p[100], x[100], y[100];
  double X,Y ;
  int num_of_points = 0;
  double terms[100];

  // must do this before you do 'almost' any other graphical tasks 
  swidth = 830 ;  sheight = 800 ;
  G_init_graphics (swidth,sheight) ;  // interactive graphics

  // clear the screen in a given color
  G_rgb (0.3, 0.3, 0.3) ; // dark gray
  G_clear () ;

  // draw a different colored rectangle where if it is clicked, 
  // the program will go forward
  G_rgb (0.0, 0.0, 1.0) ; // blue
  lowleftx = 800; lowlefty = 0 ; width = 30 ; height = 799 ;
  G_fill_rectangle (lowleftx, lowlefty, width, height) ;
  
  // Get all the points the user wants, terminate when they click above 
  // 800 (blue bar)
  G_rgb(1,1,0) ;
  for(int i = 0; i < 100; ++i){
    G_wait_click(p);
    if (p[0] > 800)
      break;
    x[i] = p[0] ; y[i] = p[1];
    G_fill_circle(x[i],y[i],2);
    ++num_of_points;
  }
 
  // LaGrange mathy stuff
  double top_product, bottom_product;
  G_rgb(1,0,0) ;
  for (X = 0; X < 800; ++X){
    Y = 0.0;
    for (int i = 0; i < num_of_points; ++i){
      // set = to 1, so we don't divide by 0
      top_product = 1.0;
      bottom_product = 1.0;
      for (int j = 0; j < num_of_points; ++j){
        // skip so we don't get a 0 in the denominator
        if (j == i){
          continue;
        }
        top_product *= (X - x[j]);
        bottom_product *= (x[i] - x[j]); 
      }
      Y += y[i] * (top_product/bottom_product);
    }
    G_point(X,Y);
  }

  int key ;   
  key =  G_wait_key() ; // pause so user can see results

  G_save_image_to_file("three_point_quadratic.xwd") ;
}
