/* Hadamard Transform
   mex function to take hadamard transform

   Usage: w = hadamard(x)
   x must be a REAL VALUED COLUMN VECTOR or MATRIX
   m = size(x,1) must be a POWER OF TWO

   Notes:
   1) This implementation uses exactly m*log2(m) additions/subtractions.
   2) This is symmetric and orthogonal. To invert, apply again and
      divide by vector length.

   Written by: Peter Stobbe, Caltech
   Email: stobbe@acm.caltech.edu
   Created: August 2008
   Edits by Stephen Becker, 2009--2014

      Note: in R2008b, Matlab added "fwht" and "ifwht" (the Fast Walsh-
          Hadamart Transform and the inverse) to its Signal Processing
          Toolbox.  With the default ordering and scaling, it's not
          equivalent to this, but you can change this with the following:
          y = length(x) * fwht( x, [], 'hadamard' );
          Then y should be the same as hadamard(x) up to roundoff.
          However, it appears that this code is faster than fwht.

 Update Stephen Becker, Feb 27 2014, fix compiling issue for Mac OS X
 Update Stephen Becker, Mar  3 2014, issue error if input data is sparse
*/

#include <stdlib.h>


/* SRB: Feb 27 2014, gcc-4.8 has problems with char16_t not being defined. 
 * This  seems to fix it
 * (and do this BEFORE including mex.h) */
/* See http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56086#c4 
 (but for, e.g., Mac w/ Xcode and Clang, this fails, so test
 for gcc. more possibilities here:
  https://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
 but clang defines GNUC too!
 http://nadeausoftware.com/articles/2012/10/c_c_tip_how_detect_compiler_name_and_version_using_compiler_predefined_macros
 */
#ifndef NO_UCHAR
#define UCHAR_OK
#endif
#if defined(__GNUC__) && !(defined(__clang__)) && defined(UCHAR_OK)
#include <uchar.h>
#endif


/* 
 y - output
 x - input
 m - length of vector
 */
void hadamard_apply_vector(double *y, double *x, unsigned m)
{
  unsigned bit, j, k;
  double temp;
   
  for (j = 0; j < m; j+=2) {
      k = j+1;
      y[j] = x[j] + x[k];
      y[k] = x[j] - x[k];
  }
  
  for (bit = 2; bit < m; bit <<= 1) {   
    for (j = 0; j < m; j++) {
        if( (bit & j) == 0 ) {
              k = j | bit;
              temp = y[j];
              y[j] = y[j] + y[k];
              y[k] = temp - y[k];
        }
    }
  }
}

void hadamard_apply_vector_inplace(double *x, unsigned m) {
   // http://www.musicdsp.org/showone.php?id=18
   unsigned log2m;
   for (log2m = 0; m > 1; m >>= 1) {
      ++log2m;
   }
   
   unsigned i,j,k;
   unsigned jbnd = 1 << log2m;
   unsigned one_shift_i, one_shift_ip1;
   unsigned jpk, jpk_one_shift_i;
   double temp;
   for (i = 0; i < log2m; ++i) {
      one_shift_i = (1 << i);
      one_shift_ip1 = (1 << (i+1));

      for (j = 0; j < jbnd; j += one_shift_ip1) {
         for (k = 0; k < one_shift_i; ++k) {
            jpk = j+k;
            jpk_one_shift_i = jpk + one_shift_i;

            temp = x[jpk];
            x[jpk] += x[jpk_one_shift_i];
            x[jpk_one_shift_i] = temp - x[jpk_one_shift_i];
         }
      }
   }
}


///* 
// y - output
// x - input
// m - length of vectors (number of rows)
// n - number of vectors (number of columns)
// */
//void hadamard_apply_matrix(double *y, double *x, unsigned m, unsigned n)
//{
//    unsigned j;
//    for(j = 0; j < n; j++) {
//        hadamard_apply_vector(y + j*m, x + j*m, m);
//    }
//}


