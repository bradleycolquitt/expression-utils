library(inline)
library(Rcpp)
sourceCpp("/media/data2/rstudio/birds/utils/subsample.cpp")

src1 = 
'
Rcpp::NumericVector dv = Rcpp::clone<Rcpp::NumericVector>(d);
'
src <- 
  '
int num = as<int>(size), x = as<int>(n);
Rcpp::NumericVector vx = Rcpp::clone<Rcpp::NumericVector>(x);
Rcpp::NumericVector pr = Rcpp::clone<Rcpp::NumericVector>(prob);
Rcpp::NumericVector rnd = rexp(x) / pr;
for(int i= 0; i<vx.size(); ++i) vx[i] = i;
std::partial_sort(vx.begin(), vx.begin() + num, vx.end(), Comp(rnd));
vx = vx[seq(0, num - 1)] + 1;
return vx;
'
incl <- 
  '
struct Comp{
  Comp(const Rcpp::NumericVector& v ) : _v(v) {}
  bool operator ()(int a, int b) { return _v[a] < _v[b]; }
  const Rcpp::NumericVector& _v;
};
'
funFast <- cxxfunction(signature(n = "Numeric", size = "integer", prob = "numeric"),
                       src, plugin = "Rcpp", include = incl)

#' faster weighted random sample
#' @references, http://stackoverflow.com/questions/15113650/faster-weighted-sampling-without-replacement
weighted_Random_Sample = function(ldata, .weights, .n) {
  key = runif(ldata) ^ (1 / .weights)
  return(which.max(key))
  #return(order(key, decreasing = T)[1:.n])
#  return(.data[order(key, decreasing=TRUE)][1:.n])
}

subsample_umi = function(d, column, target, max_val, max_count) {
  total_counts = sum(d[,column])
  if (total_counts <= target)
    return(d)
  
  d = as.data.frame(d)
  ldata = nrow(d)
  vec = as.vector(d[,column])
  for (i in total_counts:(target-1)) {
    ind = funFast(ldata, 1, vec)
    vec[ind] = vec[ind] - 1
  }
  d[,column] = vec
  return(d)
}


subsample_umi2 = function(d, column, target, max_val, max_count) {
  total_counts = sum(d[,column])
  if (total_counts <= target)
    return(d)
  
  d = as.data.frame(d)
  vec = as.vector(d[,column])
  d[,column] = subsample_umi_cpp(vec, target)
  return(d)
}

subsample_umi3 = function(d, column, target, max_val, max_count) {
 
  total_counts = sum(d[,column])
  if (total_counts <= target)
    return(d)
  
  d = as.data.frame(d)
  vec = as.vector(d[,column])
  d[,column] = subsample_umi1_cpp(vec, target)
  return(d)
}


#' Subsample dataset iteratively, Using counts in 'column' as weights
#' Intended to subsample data.frame containing UMI counts
#' @param d , data.frame
#' @param column , column to subsample
#' @param target , target total number of counts in 'column'
#'
#' @return data.frame as d with subsampled 'column'
#' @export 
#'
#' @examples
subsample_umi4 = function(d, column, target) {
  
  total_counts = sum(d[,column])
  if (total_counts <= target)
    return(d)
  
  d = as.data.frame(d)
  vec = as.vector(d[,column])
  d[,column] = subsample_umi2_cpp(vec, target, step = 10)
  return(d)
}

