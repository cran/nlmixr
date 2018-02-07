#include <RcppArmadillo.h>
#include <R_ext/Rdynload.h>
#include <Eigen/Dense>
#include "saem_class_rcpp.hpp"
#include "lin_cmt.hpp"


using namespace std;
using namespace arma;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

Function ff("sd");

RObject mat2NumMat(const mat &m) {
	RObject x = wrap( m.memptr() , m.memptr() + m.n_elem ) ;
	x.attr( "dim" ) = Dimension( m.n_rows, m.n_cols ) ;
	return x;
}

vec Ruser_function(const mat &phi_, const mat &evt_, const List &opt) {
  RObject phi, evt;
  phi = mat2NumMat(phi_);
  evt = mat2NumMat(evt_);
  NumericVector g;
  g = ff(phi, evt);
  vec yp(g);

  return yp;
}


vec user_function(const mat &phi, const mat &evt, const List &opt) {
  uvec ix;
  vec id = evt.col(0);
  mat wm;
  vec obs_time, dose_time, dose, wv;

  ix = find(evt.col(2) == 0);
  vec yp(ix.n_elem);
  double *p=yp.memptr();
  int N=id.max()+1;

  for (int i=0; i<N; i++) {
    ix = find(id == i);
    wm = evt.rows(ix);

    ix = find(wm.col(2) == 0);
    wv = wm.col(1);
    wv = wv(ix);
    const Map<MatrixXd> _obs_time(wv.memptr(), wv.n_elem, 1);
    const VectorXd obs_time(_obs_time);

    ix = find(wm.col(2) > 0);
    wv = wm.col(1);
    wv = wv(ix);
    const Map<MatrixXd> _dose_time(wv.memptr(), wv.n_elem, 1);
    const VectorXd dose_time(_dose_time);

    wv = wm.col(3);
    wv = wv(ix);
    const Map<MatrixXd> _dose(wv.memptr(), wv.n_elem, 1);
    const VectorXd dose(_dose);

    wv = wm.col(4);
    wv = wv(ix);
    const Map<MatrixXd> _Tinf(wv.memptr(), wv.n_elem, 1);
    const VectorXd Tinf(_Tinf);

  double
  CL,
  lCL,
  V,
  lV,
  KA,
  lKA;

  lCL = phi(i, 0);
  lV = phi(i, 1);
  lKA = phi(i, 2);

  CL = exp( lCL);
  V = exp( lV);
  KA = exp( lKA);
VectorXd params(4);
params(0) = CL;
params(1) = V;
params(2) = KA;
params(3) = 0;

    int no=obs_time.size();
    VectorXd g(obs_time.size());
    int ncmt=1, oral=TRUE, infusion=FALSE, parameterization=1;

	g = generic_cmt_interface(
      obs_time,
      dose_time,
      dose,
      Tinf,
      params,
      ncmt,
      oral,
      infusion,
      parameterization);

    memcpy(p, g.data(), no*sizeof(double));
    p += no;
	//cout << "ok " << i <<endl;
  }

  return yp;
}

// definition
RcppExport SEXP dopred( SEXP in_phi, SEXP in_evt, SEXP in_opt ) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< mat& >::type phi(in_phi);
    Rcpp::traits::input_parameter< mat& >::type evt(in_evt);
    List opt(in_opt);
    vec g = user_function(phi, evt, opt);
    return Rcpp::wrap(g);
END_RCPP
}

RcppExport SEXP saem_fit(SEXP xSEXP) {
BEGIN_RCPP
  List x(xSEXP);
  SAEM saem;
  saem.inits(x);

  if(x.containsElementNamed("Rfn")) {
    ff = as<Function>(x["Rfn"]);
    saem.set_fn(Ruser_function);
  } else {
    saem.set_fn(user_function);
  }

  saem.saem_fit();

  List out = List::create(
    Named("mpost_phi") = saem.get_mpost_phi(),
    Named("Gamma2_phi1") = saem.get_Gamma2_phi1(),
    Named("Plambda") = saem.get_Plambda(),
    Named("Ha") = saem.get_Ha(),
    Named("sig2") = saem.get_sig2(),
    Named("eta") = saem.get_eta(),
    Named("par_hist") = saem.get_par_hist()
  );
  out.attr("saem.cfg") = x;
  out.attr("class") = "saemFit";
  return out;
END_RCPP
}

