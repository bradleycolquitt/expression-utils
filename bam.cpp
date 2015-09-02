#include <Rcpp.h>
#include <sam.h>
using namespace Rcpp;
using namespace sam;
// [[Rcpp::depends sam]]
// [[Rcpp::export]]
void pileup(char* fname) {
  BGZF* bgzf;
  bgzf = bgzf_open(fname, 'r');
  bam_hdr_t header;
  header = bam_hdr_read(bgzf);
  std::cout << header.n_targets << std::end;
  
  
}