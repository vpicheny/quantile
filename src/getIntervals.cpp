#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int find_nth(const Rcpp::NumericVector& xa, const int middle)
{
  Rcpp::NumericVector sort(middle);
  /*std::cout << "middle=" << middle << std::endl;*/
  std::partial_sort_copy(xa.begin(), xa.end(), sort.begin(), sort.end());
  /*for (int i = 0; i < sort.size(); ++i) {*/
  /*std::cout << \' \' << i << \':\' << xa[i] << \'/\' << sort[i];*/
  /*}*/
  /*std::cout << std::endl;*/
  return std::find(xa.begin(), xa.end(), sort[sort.size() - 1]) - xa.begin() + 1;  // because R.
}

// [[Rcpp::export]]
int
c_iter_kk(const Rcpp::NumericVector& c1, const Rcpp::NumericVector& c2, const int i2, int imax, int value)
{
  int i = i2;
  for (; i < imax && c1[i] != value && c2[i] != value; ++i);
  return i + 1;
}


template <typename Scalar>
struct gen_idx {
Scalar i;
gen_idx() : i(0) {}
Scalar operator () () { return i += 1; }
};

struct comp_idx {
const Rcpp::NumericVector* ref;
comp_idx(const Rcpp::NumericVector& r) : ref(&r) {}

bool operator () (double d1, double d2) const
{
  int i1 = (int) d1;
  int i2 = (int) d2;
  return (*ref)[i1 - 1] < (*ref)[i2 - 1];
}

bool operator () (int i1, int i2) const
{
  return (*ref)[i1 - 1] < (*ref)[i2 - 1];
}
};

// [[Rcpp::export]]
Rcpp::NumericVector
c_order_nv(const Rcpp::NumericVector& vec)
{
  Rcpp::NumericVector idx(vec.size());
  std::generate(idx.begin(), idx.end(), gen_idx<double>());
  std::sort(idx.begin(), idx.end(), comp_idx(vec));
  return idx;
}

// [[Rcpp::export]]
std::vector<int>
c_order(const Rcpp::NumericVector& vec)
{
  std::vector<int> idx(vec.size());
  std::generate(idx.begin(), idx.end(), gen_idx<int>());
  std::sort(idx.begin(), idx.end(), comp_idx(vec));
  return idx;
}

// [[Rcpp::export]]
Rcpp::NumericVector
getIntervals(const Rcpp::NumericVector& vec_a, const Rcpp::NumericVector& vec_b, int k, double Imin)
{
  --k;
  // cherche argmin(-(b[k] - b) / (a[k] - a) > Imin)
  int min_idx = -1;
  double min_val = std::numeric_limits<double>::infinity();
  for (int idx = 0; idx < vec_a.size(); ++idx) {
  if (idx == k) { continue; }
  double val = (vec_b[idx] - vec_b[k]) / (vec_a[k] - vec_a[idx]);
  if (val < min_val && val > Imin) {
  min_idx = idx;
  min_val = val;
  }
  }
  return Rcpp::NumericVector::create(min_idx + 1, min_val);
}