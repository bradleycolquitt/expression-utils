#include <Rcpp.h>
using namespace Rcpp;

struct Comp{
  Comp(const Rcpp::NumericVector& v ) : _v(v) {}
  bool operator ()(int a, int b) { return _v[a] < _v[b]; }
  const Rcpp::NumericVector& _v;
};

// [[Rcpp::export]]
void fast_subsample1(IntegerVector vx, NumericVector & pr, IntegerVector & ind)
{
  Rcpp::NumericVector rnd = rexp(vx.size()) / pr;
  std::partial_sort(vx.begin(), vx.begin() + 1, vx.end(), Comp(rnd));
  ind = vx[0];
  //return ;
}

// [[Rcpp::export]]
Rcpp::NumericVector fast_subsample(int size, int n, NumericVector pr)
{
Rcpp::NumericVector vx = Rcpp::clone<Rcpp::NumericVector>(size);
Rcpp::NumericVector rnd = rexp(size) / pr;
for(int i= 0; i<vx.size(); ++i) vx[i] = i;
std::partial_sort(vx.begin(), vx.begin() + n, vx.end(), Comp(rnd));
vx = vx[seq(0, n - 1)];
return vx;
}

// [[Rcpp::export]]
void fast_subsample2(IntegerVector vx, NumericVector & pr, IntegerVector & ind, int n)
{
 // Rcpp::NumericVector vx = Rcpp::clone<Rcpp::NumericVector>(vx.size());
  Rcpp::NumericVector rnd = rexp(vx.size()) / pr;
  //for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  std::partial_sort(vx.begin(), vx.begin() + n, vx.end(), Comp(rnd));
  //ind = vx[seq(0,n)];
  for (int j=0; j<n; j++) {
    pr[vx[j]]--;
  }
 // vx = vx[seq(0, n - 1)];
  //return vx;
}
// [[Rcpp::export]]
Rcpp::NumericVector subsample_umi_cpp(Rcpp::NumericVector vec, int target) 
  {
  int total_counts = sum(vec);
  //Rcout << total_counts << std::endl;
  int ind = 0;
  Rcpp::NumericVector pr = Rcpp::NumericVector(vec.size());
  Rcpp::IntegerVector vx = Rcpp::clone<Rcpp::IntegerVector>(vec.size());
  for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  for (int i=total_counts; i>target; i--) 
    {
//pr =  vec / i;
// pr[pr<0] = 0;
    //Rcout << pr << std::endl;
    ind = fast_subsample(vec.size(), 1, vec)[0];
    vec[ind]--;
    //Rcout << vec << std::endl;
  }
  return vec;
}

// [[Rcpp::export]]
Rcpp::NumericVector subsample_umi1_cpp(Rcpp::NumericVector vec, int target) 
{
  int total_counts = sum(vec);
  //Rcout << total_counts << std::endl;
  //int ind = 0;
  Rcpp::NumericVector pr = Rcpp::NumericVector(vec.size());
  Rcpp::IntegerVector vx = Rcpp::clone<Rcpp::IntegerVector>(vec.size());
  //Rcpp::IntegerVector vx_tosort = Rcpp::clone<Rcpp::IntegerVector>(vec.size());
  for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  
  Rcpp::IntegerVector ind = Rcpp::IntegerVector(1);
  for (int i=total_counts; i>target; i--) 
  {
    //pr =  vec / i;
    // pr[pr<0] = 0;
    //Rcout << pr << std::endl;
    fast_subsample1(vx, vec, ind);
    vec[ind[0]]--;
    //Rcout << vec << std::endl;
  }
  return vec;
}

// [[Rcpp::export]]
Rcpp::NumericVector subsample_umi2_cpp(Rcpp::NumericVector vec, int target, int step) 
{
  int total_counts = sum(vec);
  //Rcout << total_counts << std::endl;
  //int ind = 0;
  //Rcpp::NumericVector vec_local = Rcpp::clone<Rcpp::NumericVector>(vec.size());
  Rcpp::NumericVector pr = Rcpp::clone<Rcpp::NumericVector>(vec);
  Rcpp::IntegerVector vx = Rcpp::clone<Rcpp::IntegerVector>(vec.size());
  //Rcpp::IntegerVector vx_tosort = Rcpp::clone<Rcpp::IntegerVector>(vec.size());
  for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  
  Rcpp::IntegerVector ind = Rcpp::IntegerVector(step);
  for (int i=total_counts; i>target; i=i-step) 
  {
    //pr =  vec / i;
    // pr[pr<0] = 0;
    //Rcout << pr << std::endl;
    fast_subsample2(vx, pr, ind, step);
    //vec_local[]--;
    //vec_local[ind[seq(0,step)]] = vec_local[ind[seq(0,step)]] - 1;
    //Rcout << vec << std::endl;
  }
  return pr;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
fast_subsample(42, 1, rep(1, times=42)/42)
*/

/*** R
d = 1:10
subsample_umi_cpp(d, 10)
*/

/*** R
d = 1:10
subsample_umi2_cpp(d, 10, 2)
*/

