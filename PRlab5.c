#include "FPToolkit.c"
#include <stdbool.h>

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

int main()
{
  int swidth, sheight;
  int degree;
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
  G_rgb(1,0,0) ;
  
  // intialize
  double matrix[M][M+1];
  for (int i = 0; i < M; ++i)
    for (int j = 0; j <= M; ++j)
      matrix[i][j] = 0;
  int num_of_vars = (num_of_points - 1) * 2;
  int current = 0;
  bool first = true;
  // 6 points will yield 5 splines (10 variables)
  // A and B are the line variables
  // special cases
  matrix[0][0] = 2;
  matrix[0][1] = -2 * (x[current - 1] - x[current]);
  matrix[0][num_of_vars] = 0;

  for (int i = 1; i < num_of_vars; ++i){
    for (int j = 0; j < num_of_vars + 1; ++j){
      if (first){
        
      }
      matrix[i][j] = 0;
    }

  }

  print_matrix(matrix, num_of_vars);

  // graph the function using Y = ax^0 + bx^1 + cx^2...
  /*
  for (X = 0; X < 800; ++X){
    for (int i = 0; i < degree +1; ++i){
      // get the Y value for the best fit
    }
    G_point(X,Y);
  }
  int key ;   
  key =  G_wait_key() ; // pause so user can see results

  G_save_image_to_file("Best_fit_poly.xwd") ;
  */
}
