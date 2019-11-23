
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

#define M 20

// Using david's gauss_elim, just in case
void print_matrix (double m[M][M+1], int n)
{
  int r,c ;

  for (r = 0 ; r < n ; r++) {
    for (c = 0 ; c <= n ; c++) {
      printf("%.16lf ",m[r][c]) ;
    }
    printf("\n") ;

  }
  printf("\n") ;
}

int gaussian_elimination(double m[M][M+1], int n, double x[M])
{
  int i,j,k ;
  double v,sum ;
  int max_row,c ;

  // reduce matrix to upper triangular form 

  for (j = 0 ; j < n-1 ; j++) {

    //print_matrix(m,n) ;

    /*
    // search from row j+1 down to find the largest magnitude
    max_row = j ;
    for (k = j+1 ; k < n ; k++) {
    if (fabs(m[k][j]) > fabs(m[max_row][j])) { max_row = k ; }
    }
    if (max_row != j) {
    // swap rows
    for(c = j ; c <= n ; c++) {
    v = m[j][c] ; m[j][c] = m[max_row][c] ; m[max_row][c] = v ;
    }
    }
    */

    if (m[j][j] == 0) return 0 ; 
    // this now means there does not exist a unique soln
    // have to think about this ...
    // consider the subproblem underneath this
    // it has either a unique soln, is inconsistent, or is underdetermined
    // Now we are adding a new equation with 0*new_variable + ...old vars...
    // Since 0 is multiplying the new_variable, we can view this
    // as an extra equation being added to the subproblem which is
    // already full of equations.  
    // case 1: if the subproblem had a unique solution then this solution
    //         either works or not with this added equation.  if it does
    //         not work, then we have an inconsistent system
    //         if it DOES work, then throwing the extra dimension in
    //         creates an UNDERDETERMINED situation as any value of
    //         xnew will still make things work
    // case 2: the subproblem is inconsistent. so no set of values
    //         of xsub can make the sub problem work and that will
    //         still be true, so new system is also inconsistenet
    // case 3: the subproblem is underdetermined.  throwing in the
    //         new equation could create a unique solution in the
    //         subspace but still leaves any choice for xnew so 
    //         the larger system is still underdetermined



    for (k = j+1 ; k < n ; k++) {

      v = m[k][j] / m[j][j] ;
      for (i = 0 ; i <= n ; i++) {
        m[k][i] = m[k][i] - v*m[j][i] ;
      }

    }


  }

  // Output the upper triangular form
  //print_matrix(m,n) ;


  //////////////////////////////////////////////////////

  // Now do the back substitution 
  for (j = n - 1 ; j >= 0 ; j--) {

    sum = 0.0 ;
    for (k = j + 1 ; k < n ; k++) {
      sum += (m[j][k] * x[k]) ;
    }

    if (m[j][j] == 0) return 0 ;

    x[j] = (m[j][n] - sum)/m[j][j] ;
  }

  return 1 ;
}


int main()
{
  int swidth, sheight;
  int degree;
  double lowleftx, lowlefty, width, height ;
  double p[100], x[100], y[100];
  double X,Y ;
  int num_of_points = 0;
  double terms[100];
  // use gline() maybe?
  double oldX, oldY;
     printf("What degree polynomial would you like? ");
     scanf("%d", &degree);
     printf("Degree: %d\n", degree);
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
  G_rgb(1,0,0) ;
  /*
  int num_of_points = 4;
  degree = 2;
  // (3,4) 
  x[0] = 3;
  y[0] = 4;
  // (4,6)
  x[1] = 4;
  y[1] = 6;
  // (6,5)
  x[2] = 6;
  y[2] = 5;
  // (7,8);
  x[3] = 7;
  y[3] = 8;
*/
  // Polynomial of best fit
  double matrix[M][M+1];
  double coefficients[degree+1];
  int first_term_power = 0;
  double sum = 0.0;
  matrix[degree][degree+1] = 999;
  for (int i = 0; i <= degree; ++i){
    for (int j = 0; j <= degree +1; ++j){
      if (j == degree + 1){
        for (int k = 0; k < num_of_points; ++k){
          if (i == 0){
            sum += y[k];
            //printf("sum: %lf\n", sum);
          }
          else
            sum += y[k] * pow(x[k], i);
          //printf("(%lf,%lf) ", x[k], y[k]);
        }
        matrix[i][j] = sum;
        //printf("matrix[%i][%i] = %lf\n", i, j, sum);
        sum = 0.0;
      }
      else{
        for (int k = 0; k < num_of_points; ++k){
          sum += pow(x[k], j+i);
        }
        matrix[i][j] = sum;
        //printf("matrix[%i][%i] = %lf\n", i, j, sum);
        sum = 0.0;
      }
    }
  }
  print_matrix (matrix, degree +1);
  int success = gaussian_elimination(matrix, degree+1, coefficients);
  printf("\nCoefficients array:\n");
  for(int i = 0; i <= degree; ++i){
    printf("%lf ", coefficients[i]);
  }
  printf("\n");


  for (X = 0; X < 800; ++X){
    Y = 0.0;
    for (int i = 0; i < degree +1; ++i){
      Y += pow(X, i) * coefficients[i];
    }
    /*
    oldX = X;
    oldY = Y;
    if (X == 0)
      continue;
    G_line(oldX, oldY, X, Y);
    */
    G_point(X,Y);
  }
  int key ;   
  key =  G_wait_key() ; // pause so user can see results

  G_save_image_to_file("Best_fit_poly.xwd") ;
}
