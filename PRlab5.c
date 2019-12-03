#include "FPToolkit.c"
#include <stdbool.h>

#define M 20

// Using david's gauss_elim, just in case
void print_matrix (double m[M][M+1], int n)
{
  int r,c ;

  for (r = 0 ; r < n ; r++) {
    for (c = 0 ; c <= n ; c++) {
      printf("%.2f ",m[r][c]) ;
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
  double coefficients[num_of_vars];
  int current = 1;
  int start = 0;
  bool first = true;
  // 6 points will yield 5 splines (10 variables) so 11 columns
  // A and B are the line variables

  matrix[0][0] = 2;
  matrix[0][1] = -2 * (x[1] - x[0]);
  matrix[0][num_of_vars] = 0;

  for (int i = 1; i < num_of_vars - 1; ++i, ++start){
    for (int j = start; j < num_of_vars + 1; ++j){
      if (first){
        if (j == start){
          matrix[i][j] = x[current] - x[current-1];
        }
        else if (j == start + 1){
          matrix[i][j] = pow(x[current] - x[current-1], 2);
        }
        else if (j == start + 2){
          matrix[i][j] = x[current+1] - x[current];
        }
        else if (j == num_of_vars){
          matrix[i][j] = ((y[current+1]-y[current])/(x[current+1]-x[current])) - ((y[current]-y[current-1])/(x[current]-x[current-1]));
        }
        else{
          continue;
        }
      }
      else{
        if (j == start){
          matrix[i][j] = pow(x[current] - x[current-1], 2);
        }
        else if (j == start + 1){
          matrix[i][j] = -(x[current + 1] - x[current -1]);
        }
        else if (j == start + 2){
          matrix[i][j] = (x[current] - x[current -1]) * (x[current +1] - x[current]);
        }
        else if (j == num_of_vars){
          matrix[i][j] = -matrix[i-1][j];
        }
        else 
          continue;
      }
    }
    if (first == true)
      first = false;
    else {
      first = true;
      ++current;
    }
  }

  matrix[num_of_vars-1][num_of_vars-2] = 1;
  matrix[num_of_vars-1][num_of_vars-1] = -2 * (x[num_of_points-1] - x[num_of_points-2]);
  matrix[num_of_vars-1][num_of_vars] = 0;
  print_matrix(matrix, num_of_vars);
  // coefficients now stored in matrix 'coefficients'
  int success = gaussian_elimination(matrix, num_of_vars, coefficients);
  if (success == 0){
    printf("Gaussian failed.\n");
    return 0;
  }

  for (int i = 0; i < num_of_points; ++i)
    printf("x[%i]: %lf\n", i, x[i]);

  current = 0;
  int coeff_current = 0;
  for (X = 0; X < 800; ++X){
    if(X > x[0]){
      if (X > x[num_of_points-1])
        break;
      Y = y[current] + 
        ((y[current+1]-y[current])/(x[current+1]-x[current]))*(X-x[current]) + 
        (coefficients[coeff_current]*(X-x[current])*(X-x[current+1])) + 
        (coefficients[coeff_current+1]*pow(X-x[current],2)*(X-x[current+1]));
      G_point(X,Y);
      if(X == x[current+1]){
        ++current;
        coeff_current+=2;
      }
    }
  }
  int key ;   
  key =  G_wait_key() ; // pause so user can see results
}
