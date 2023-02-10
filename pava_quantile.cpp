#include <Rcpp.h>
#include <vector>           // for vector class
#include <algorithm>        // for sort function
#include <numeric>          // for iota function

using namespace Rcpp;
using namespace std;

// for each group we keep y values sorted to calculate quantiles in O(1)
// Is it faster to search linearly or binary search tree
// -> for small groups linear is more efficient, following
// (https://stackoverflow.com/questions/36041553/sorting-one-element-in-a-sorted-container)

// it is not clear how many group we need: should g_val and g_i_start be initialized
// statically with n elements or should it be able to grow dynamically?

// use c++ convention: lower bound inclusive, upper bound exclusive
// i++ means i is incremented after current line is executed (postfix)

// [[Rcpp::export]]
NumericMatrix pav_quantile_cpp(NumericVector x, NumericVector y, NumericVector alpha) {
  if (x.size() != y.size() || alpha.size() != 1) {
    return -1;
  }
  int n = x.size();
  double a = alpha[0];
  vector<int> ord(n);               // order of x values
  vector<double> y_arr(n);
  vector<int> g_i_start(n);         // start index of each group
  vector<double> g_val(n);          // quantile of each group
  vector<double> sorted(n);         // is used temporally for sorting
  int g_curr = 0;
  int i, k, k1, k2, size;
  vector<double>::iterator i_insert;

  // determine order of x
  iota(ord.begin(), ord.end(), 0);                   // fill vector with 1:n
  sort(ord.begin(), ord.end(),      // sort with x and in case of ties with y
       [&](int A, int B) -> bool { return (x[A] < x[B] || (x[A] == x[B] && y[A] > y[B])); });

  // now walk through y-values and pool adjacent violators
  for (i = 0; i < n; i++) {
    // create new group plus pool violations looking forward -------------------
    g_i_start[g_curr] = i;
    y_arr[i] = y[ord[i]];
    g_val[g_curr] = y_arr[i];
    while (i < n - 1 && g_val[g_curr] > y[ord[i + 1]]) {
      // sort y[ord[i+1]] into current group
      i++;
      i_insert = lower_bound(y_arr.begin() + g_i_start[g_curr], y_arr.begin() + i, y[ord[i]]);
      // i_insert points where new value has to be sorted in: copy succeeding values one back
      copy_backward(i_insert, y_arr.begin() + i, y_arr.begin() + i + 1);
      *i_insert = y[ord[i]];
      g_val[g_curr] = y_arr[g_i_start[g_curr] + (int) (a * (i + 1 - g_i_start[g_curr]))];
    }

    // pool violations looking backward ----------------------------------------
    while (g_curr > 0 && g_val[g_curr] < g_val[g_curr - 1]) {
      // [g_i_start[g_curr-1], g_i_start[g_curr]) and [g_i_start[g_curr], i) are sorted
      // --> merge them by traversing both sections simultaneously
      size = i + 1 - g_i_start[g_curr - 1];  // cannot do it inplace
      k1 = g_i_start[g_curr - 1];
      k2 = g_i_start[g_curr];
      for (k = 0; k < size; k++) {
		// all elements from 2nd list used or 1st list has smaller element -> put it next
        if (k2 > i || (k1 < g_i_start[g_curr] && y_arr[k1] < y_arr[k2])) {
          sorted[k] = y_arr[k1++];
        } else {
          sorted[k] = y_arr[k2++];
        }
      }
	  g_i_start[g_curr] = 0;		// delete old group
      g_curr--; // we merge two groups, so there is one less
      for (k = 0; k < size; k++) {
        y_arr[g_i_start[g_curr] + k] = sorted[k];          // copy sorted array
      }
      // determine quantile (casting to int corresponds flooring function)
      g_val[g_curr] = y_arr[g_i_start[g_curr] + (int)(size * a)];
    }
    g_curr++;
  }
  // build output together -----------------------------------------------------
  NumericMatrix out(n, 2);
  int i_next = 0;
  for (i = 0; i < n; i++) {
    out(i, 0) = x[ord[i]];
    out(i, 1) = g_val[i_next];
    // if we are in next group switch value
    if (i_next < n - 1 && g_i_start[i_next + 1] > 0 && i + 1 >= g_i_start[i_next + 1]) {
      i_next++;
    }
  }
  return out;
}
