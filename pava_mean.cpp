#include <Rcpp.h>
#include <vector>           // for vector class
#include <algorithm>        // for sort function
#include <numeric>          // for iota function

using namespace Rcpp;
using namespace std;

// important: coinciding x-values have to be resolved immediately
// made me some trouble finding that error

inline void pool_with_prev(int i, vector<double> &g_val, vector<int> &g_count) {
    int new_count = g_count[i] + g_count[i - 1];
    // update mean and count
    g_val[i-1] = (g_count[i] * g_val[i] + g_count[i-1] * g_val[i-1]) / new_count;
    g_count[i-1] = new_count;
}

// [[Rcpp::export]]
NumericMatrix pav_mean_cpp(NumericVector x, NumericVector y) {
  if (x.size() != y.size()) {
    return -1;
  }
  int n = x.size();
  int i;
  vector<int> order(n);         // order of x values
  vector<double> x_sorted(n);
  vector<int> g_count(n);       // size of each group
  vector<double> g_val(n);         // value of each group

  // determine order of x
  iota(order.begin(), order.end(), 0);                   // fill vector with 1:n
  sort(order.begin(), order.end(),
       [&](int A, int B) -> bool { return x[A] < x[B]; });  // lambda expression
  for (i = 0; i < n; i++) {
    x_sorted[i] = x[order[i]];
  }

  // now walk through y-values and pool adjacent violators
  int g_curr = 0;
  for (i = 0; i < n; i++) {
    // create new group
    g_val[g_curr] = y[order[i]];
    g_count[g_curr] = 1;
    while (i < n - 1 && x_sorted[i] == x_sorted[i + 1]) { // coinciding x-values?
        i++;
        g_val[g_curr] += y[order[i]];
        g_count[g_curr]++;
    }
    g_val[g_curr] /= g_count[g_curr];                     // average

    if (g_curr > 0 && g_val[g_curr] < g_val[g_curr - 1]) { // violations?
      pool_with_prev(g_curr, g_val, g_count);
      // new violations before? --> pool them!
      while ((g_curr > 1) && (g_val[g_curr - 1] < g_val[g_curr - 2])) {
	    g_curr--;  	// decrease beforehand as we have checked for j-1 and j-2
        pool_with_prev(g_curr, g_val, g_count);
      }
    } else {                    // otherwise just move on
      g_curr++;
    }
  }
  // build output together
  NumericMatrix out(n, 2);
  int i_next = 0;
  for (i = 0; i < n; i++) {
    out(i, 0) = x_sorted[i];
    out(i, 1) = g_val[i_next];
    g_count[i_next]--;              // reduce count
    if (g_count[i_next] == 0) {     // if count is zero, switch to next value
      i_next++;
    }
  }
  return out;
}
