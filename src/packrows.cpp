
#include "Rcpp.h"
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::IntegerVector
packrows_cpp(IntegerVector tmin, IntegerVector tmax) {
  int n = tmin.size();
  int nr = n;  // counter for completed genes
  IntegerVector trow(n);
  int j = 1;
  
  while(nr > 0)
  {
    int jrow = 0;
    IntegerVector jset(n);  // accumulate placements in each row
    
    for(int i = 0; i < n; ++i) {
      if (trow[i] == 0) {  // which need placement
        if (jrow == 0) {  // place 1st
          trow[i] = j;
          jset[0] = i;
          jrow = 1;
          nr--;
        } else {
          // check existing placements
          for (int k = 0; k < jrow; ++k) {
            if (tmin[i] < tmax[jset[k]] && tmax[i] > tmin[jset[k]]) {
              break;
            }
            if (k == jrow - 1) {
              // no clashes
              trow[i] = j;
              jset[jrow++] = i;
              nr--;
            }
          }
        }
      }
    }
    j++;
  }
  
  return trow;
}
