//[[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// ////////////////////////////////////////////////////////////////////////

//' Great circle distance matrix.
//'
//' Given a set of longitude/latitude locations, \code{rdist_earth_cpp}
//' computes the Great circle (geographic) distance matrix among all pairings.
//'
//' Surprisingly this is all done efficiently in R by dot products of the
//' direction cosines.  This function is a (partial!) port to C++ of D. Nychka's
//' \code{rdist.earth} function from the fields R package.
//'
//' @param x Matrix of first set of lon/lat coordinates; first column is
//' the longitudes and the second is the latitudes.
//' @params miles If true, distances are in statute miles. If false, distances
//' are in kilometers.
//' @params R Radius to use for sphere to find spherical distances. If NULL,
//' the radius is either in miles or km depending on the values of the miles
//' argument. If R = 1, then distances are of course in radians.
// [[Rcpp::export]]
arma::mat rdist_earth_cpp(NumericMatrix x1, NumericMatrix x2,
                          bool miles = true){

  double R;
  (miles) ? ( R =  3963.34 ) : (R  = 6378.388) ;
  const double pi = 3.14159265358979323846264;
  int N1 = x1.nrow(); //std::cout << N1 << '\n';
  int N2 = x2.nrow(); //std::cout << N2 << '\n';

  NumericVector coslat1(N1), sinlat1(N1), coslon1(N1), sinlon1(N1);
  NumericVector coslat2(N2), sinlat2(N2), coslon2(N2), sinlon2(N2);
  arma::mat M1(N1,3), M2(N2,3);

  for(int i = 0; i < N1; i++){
    coslat1(i) = cos((x1(i, 1) * pi)/180);
    sinlat1(i) = sin((x1(i, 1) * pi)/180);
    coslon1(i) = cos((x1(i, 0) * pi)/180);
    sinlon1(i) = sin((x1(i, 0) * pi)/180);

    M1(i,0) = coslat1(i) * coslon1(i);
    M1(i,1) = coslat1(i) * sinlon1(i);
    M1(i,2) = sinlat1(i);
  }

  for(int i = 0; i < N2; i++){
    coslat2(i) = cos((x2(i, 1) * pi)/180);
    sinlat2(i) = sin((x2(i, 1) * pi)/180);
    coslon2(i) = cos((x2(i, 0) * pi)/180);
    sinlon2(i) = sin((x2(i, 0) * pi)/180);

    M2(i,0) = coslat2(i) * coslon2(i);
    M2(i,1) = coslat2(i) * sinlon2(i);
    M2(i,2) = sinlat2(i);
  }

  arma::mat pp = M1 * M2.t();

  arma::mat D = R* acos(pp);
  D.elem(find(abs(pp)>1)) =
    R * acos(1 * sign(pp.elem(find(abs(pp)>1))));

  return(D);
}

