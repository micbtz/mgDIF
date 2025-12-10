// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;



// [[Rcpp::export]]
double ProfLik_unc_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                        arma::mat aj1T, arma::mat bj1T, arma::mat var_aj1T, 
                        arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec alpha;
  arma::vec beta;
  arma::vec A = arma::ones(numforms);
  arma::vec B = arma::zeros(numforms);
  arma::vec Ashort = par.subvec(0, numforms - 2);
  arma::vec Bshort = par.subvec(numforms - 1, 2 * numforms - 3);
  A.elem(notbase) = Ashort; 
  B.elem(notbase) = Bshort; 
  arma::mat aj1T_res = aj1T;
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    alpha = par.subvec(2 * numforms - 2, 2 * numforms - 3 + numvar);
    beta = par.subvec(2 * numforms + numvar - 2, 2 * numforms + 2 * numvar - 3);
    for (i = 0; i < numvar; i ++)
    {
      arma::mat Xmat_i = Xmat(i);
      aj1T_res -= alpha[i] * Xmat_i;
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Amat = repmat(A.t(), aj1T.n_rows, 1);
  arma::mat numas = aj1T_res % Amat / var_aj1T;
  arma::mat denas = pow(Amat, 2) / var_aj1T;
  numas.replace(arma::datum::nan, 0);
  denas.replace(arma::datum::nan, 0);
  arma::vec numas1 = sum(numas, 1);
  arma::vec denas1 = sum(denas, 1);
  arma::vec as = numas1 / denas1;
  arma::mat Amat_as = as * A.t();
  arma::mat lla_jt = pow(aj1T_res - Amat_as, 2) / var_aj1T;
  double lla = accu(lla_jt.elem(arma::find_finite(lla_jt)));
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res % Amat + Bmat) / (var_bj1T % pow(Amat, 2));
  arma::mat denbs = 1 / (pow(Amat, 2) % var_bj1T);
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  arma::mat Bmat_bs_Amat = Bmat_bs / Amat;
  arma::mat llb_jt = pow(bj1T_res + Bmat_bs_Amat, 2) / var_bj1T;
  double llb = accu(llb_jt.elem(arma::find_finite(llb_jt)));
  double ll = lla + llb;
  return ll;
}


// [[Rcpp::export]]
double ProfLik_unc_1PL_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                            arma::mat bj1T, arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec beta;
  arma::vec B = arma::zeros(numforms);
  arma::vec Bshort = par.subvec(0, numforms - 2);
  B.elem(notbase) = Bshort; 
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    beta = par.subvec(numforms - 1, numforms + numvar - 2);
    for (i = 0; i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res + Bmat) / var_bj1T;
  arma::mat denbs = 1 / var_bj1T;
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  arma::mat llb_jt = pow(bj1T_res + Bmat_bs, 2) / var_bj1T;
  double llb = accu(llb_jt.elem(arma::find_finite(llb_jt)));
  return llb;
}


// [[Rcpp::export]]
List res_ProfLik_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                      arma::mat aj1T, arma::mat bj1T, arma::mat var_aj1T, 
                      arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec alpha;
  arma::vec beta;
  arma::vec A = arma::ones(numforms);
  arma::vec B = arma::zeros(numforms);
  arma::vec Ashort = par.subvec(0, numforms - 2);
  arma::vec Bshort = par.subvec(numforms - 1, 2 * numforms - 3);
  A.elem(notbase) = Ashort; 
  B.elem(notbase) = Bshort; 
  arma::mat aj1T_res = aj1T;
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    alpha = par.subvec(2 * numforms - 2, 2 * numforms - 3 + numvar);
    beta = par.subvec(2 * numforms + numvar - 2, 2 * numforms + 2 * numvar - 3);
    for (i = 0; i < numvar; i ++)
    {
      arma::mat Xmat_i = Xmat(i);
      aj1T_res -= alpha[i] * Xmat_i;
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Amat = repmat(A.t(), aj1T.n_rows, 1);
  arma::mat numas = aj1T_res % Amat / var_aj1T;
  arma::mat denas = pow(Amat, 2) / var_aj1T;
  numas.replace(arma::datum::nan, 0);
  denas.replace(arma::datum::nan, 0);
  arma::vec numas1 = sum(numas, 1);
  arma::vec denas1 = sum(denas, 1);
  arma::vec as = numas1/denas1;
  arma::mat Amat_as = as * A.t();
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res % Amat + Bmat) / (var_bj1T % pow(Amat, 2));
  arma::mat denbs = 1 / (pow(Amat, 2) % var_bj1T);
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  arma::mat Bmat_bs_Amat = Bmat_bs / Amat;
  
  arma::mat res_a = aj1T_res - Amat_as;
  arma::mat res_b = bj1T_res + Bmat_bs_Amat;
  Rcpp::List out = List::create(Named("res_a") = res_a, Named("res_b") = res_b);
  return out;
}

// [[Rcpp::export]]
List res_ProfLik_1PL_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                          arma::mat bj1T, arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec beta;
  arma::vec B = arma::zeros(numforms);
  arma::vec Bshort = par.subvec(0, numforms - 2);
  B.elem(notbase) = Bshort; 
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    beta = par.subvec(numforms - 1, numforms + numvar - 2);
    for (i = 0; i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res + Bmat) / var_bj1T;
  arma::mat denbs = 1 / var_bj1T;
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  arma::mat res_b = bj1T_res + Bmat_bs;
  Rcpp::List out = List::create(Named("res_b") = res_b);
  return out;
}


// [[Rcpp::export]]
List calc_asbs_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                    arma::mat aj1T, arma::mat bj1T, arma::mat var_aj1T, 
                    arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec alpha;
  arma::vec beta;
  arma::vec A = arma::ones(numforms);
  arma::vec B = arma::zeros(numforms);
  arma::vec Ashort = par.subvec(0, numforms - 2);
  arma::vec Bshort = par.subvec(numforms - 1, 2 * numforms - 3);
  A.elem(notbase) = Ashort; 
  B.elem(notbase) = Bshort; 
  arma::mat aj1T_res = aj1T;
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    alpha = par.subvec(2 * numforms - 2, 2 * numforms - 3 + numvar);
    beta = par.subvec(2 * numforms + numvar - 2, 2 * numforms + 2 * numvar - 3);
    for (i = 0;i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      aj1T_res -= alpha[i] * Xmat_i;
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Amat = repmat(A.t(), aj1T.n_rows, 1);
  arma::mat numas = aj1T_res % Amat / var_aj1T;
  arma::mat denas = pow(Amat, 2) / var_aj1T;
  numas.replace(arma::datum::nan, 0);
  denas.replace(arma::datum::nan, 0);
  arma::vec numas1 = sum(numas, 1);
  arma::vec denas1 = sum(denas, 1);
  arma::vec as = numas1/denas1;
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res % Amat + Bmat) / (var_bj1T % pow(Amat, 2));
  arma::mat denbs = 1 / (pow(Amat, 2) % var_bj1T);
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  Rcpp::List out = List::create(Named("as") = as, Named("bs") = bs);
  return out;
}


// [[Rcpp::export]]
List calc_bs_1PL_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                      arma::mat bj1T, arma::mat var_bj1T, List Xmat)
{
  int i;
  arma::vec beta;
  arma::vec B = arma::zeros(numforms);
  arma::vec Bshort = par.subvec(0, numforms - 2);
  B.elem(notbase) = Bshort; 
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    beta = par.subvec(numforms - 1, numforms + numvar - 2);
    for (i = 0;i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res + Bmat) / var_bj1T;
  arma::mat denbs = 1 / var_bj1T;
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  Rcpp::List out = List::create(Named("bs") = bs);
  return out;
}


// [[Rcpp::export]]
arma::vec derProfLik_unc_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                              arma::mat aj1T, arma::mat bj1T, arma::mat var_aj1T, 
                              arma::mat var_bj1T, List Xmat)
{
  int i, s;
  arma::vec alpha;
  arma::vec beta;
  arma::vec A = arma::ones(numforms);
  arma::vec B = arma::zeros(numforms);
  arma::vec Ashort = par.subvec(0, numforms - 2);
  arma::vec Bshort = par.subvec(numforms - 1, 2 * numforms - 3);
  A.elem(notbase) = Ashort; 
  B.elem(notbase) = Bshort; 
  arma::mat aj1T_res = aj1T;
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    alpha = par.subvec(2 * numforms - 2, 2 * numforms - 3 + numvar);
    beta = par.subvec(2 * numforms + numvar - 2, 2 * numforms + 2 * numvar - 3);
    for (i = 0;i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      aj1T_res -= alpha[i] * Xmat_i;
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Amat = repmat(A.t(), aj1T.n_rows, 1);
  arma::mat numas = aj1T_res % Amat / var_aj1T;
  arma::mat denas = pow(Amat, 2) / var_aj1T;
  numas.replace(arma::datum::nan, 0);
  denas.replace(arma::datum::nan, 0);
  arma::vec numas1 = sum(numas, 1);
  arma::vec denas1 = sum(denas, 1);
  arma::vec as = numas1/denas1;
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res % Amat + Bmat) / (var_bj1T % pow(Amat, 2));
  arma::mat denbs = 1 / (pow(Amat, 2) % var_bj1T);
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Amat_as = as * A.t();
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  arma::mat Bmat_bs_Amat = Bmat_bs / Amat;
  
  arma::mat num_der_as_A = (aj1T_res - 2 * Amat_as) / var_aj1T;
  arma::mat der_as_A = num_der_as_A;
  der_as_A.each_col() /= denas1;
  
  arma::mat num_der_bs_A = ( - bj1T_res - 2 * Bmat_bs_Amat) / (var_bj1T % pow(Amat, 2));
  arma::mat der_bs_A = num_der_bs_A;
  der_bs_A.each_col() /= denbs1;
  arma::mat num_der_bs_B = 1 / (var_bj1T % pow(Amat, 2));
  arma::mat der_bs_B = num_der_bs_B;
  der_bs_B.each_col() /= denbs1;
  arma::vec der_ll = arma::zeros(2 * numforms + 2 * numvar - 2);
  
  int count = 1;
  for (s = 0; s < numforms; s++ )
  {
    if (any(notbase == s))
    {
      arma::mat der_as_A_A = repmat(der_as_A.col(s), 1, numforms);
      der_as_A_A %= Amat;
      arma::mat der_bs_A_A = repmat(der_bs_A.col(s), 1, numforms);
      der_bs_A_A /= Amat;
      arma::mat der_bs_B_A = repmat(der_bs_B.col(s), 1, numforms);
      der_bs_B_A /= Amat;
      
      arma::mat der_lla_jt_As =  - 2 * ((aj1T_res) - Amat_as) % der_as_A_A / var_aj1T;
      arma::mat der_llb_jt_As =  - 2 * ((bj1T_res) + Bmat_bs_Amat )/var_bj1T % der_bs_A_A;
      
      double der_ll_As = accu(der_lla_jt_As.elem(arma::find_finite(der_lla_jt_As))) + 
        accu(der_llb_jt_As.elem(arma::find_finite(der_llb_jt_As)));
      
      arma::vec tobeadded =  - (aj1T_res.col(s) - A[s] * as) % as/var_aj1T.col(s) + 
        (bj1T_res.col(s) + Bmat_bs_Amat.col(s))/var_bj1T.col(s) % ( - Bmat_bs.col(s)/pow(A[s], 2));
      
      der_ll_As += 2 * accu(tobeadded.elem(arma::find_finite(tobeadded)));
      der_ll[count - 1] = der_ll_As;
      
      arma::mat der_llb_jt_Bs =  - 2 * ((bj1T_res) + Bmat_bs_Amat)/var_bj1T % der_bs_B_A;
      der_llb_jt_Bs.replace(arma::datum::nan, 0);
      double der_ll_Bs = accu(der_llb_jt_Bs);
      arma::vec tobeadded1 = (bj1T_res.col(s) - (bs - B[s])/A[s])/(var_bj1T.col(s) * A[s]);
      der_ll_Bs += 2 * accu(tobeadded1.elem(arma::find_finite(tobeadded1)));
      der_ll[count + numforms - 2] = der_ll_Bs;
      count += 1;
    }
  }
  
  arma::mat resa_jt = (aj1T_res - Amat_as) / var_aj1T;
  arma::mat resb_jt = (bj1T_res + Bmat_bs_Amat) / var_bj1T;
  
  if (numvar > 0)
  {
    for (i = 0; i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      arma::mat resaX =  - 2 * resa_jt % Xmat_i;
      arma::mat resbX =  - 2 * resb_jt % Xmat_i;
      der_ll[2 * numforms - 2 + i] = accu(resaX.elem(arma::find_finite(resaX)));
      der_ll[2 * numforms - 2 + i + numvar] = accu(resbX.elem(arma::find_finite(resbX)));
    }
  }
  return der_ll;
}



// [[Rcpp::export]]
arma::vec derProfLik_unc_1PL_Rcpp(arma::vec par, int numforms, arma::uvec notbase, 
                                  arma::mat bj1T, arma::mat var_bj1T, List Xmat)
{
  int i, s;
  arma::vec beta;
  arma::vec B = arma::zeros(numforms);
  arma::vec Bshort = par.subvec(0, numforms - 2);
  B.elem(notbase) = Bshort; 
  arma::mat bj1T_res = bj1T;
  int numvar = Xmat.size();
  if (numvar > 0)
  {
    beta = par.subvec(numforms - 1, numforms + numvar - 2);
    for (i = 0; i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      bj1T_res -= beta[i] * Xmat_i;
    }
  }
  arma::mat Bmat = repmat(B.t(), bj1T.n_rows, 1);
  arma::mat numbs = (bj1T_res + Bmat) / var_bj1T;
  arma::mat denbs = 1 / var_bj1T;
  numbs.replace(arma::datum::nan, 0);
  denbs.replace(arma::datum::nan, 0);
  arma::vec numbs1 = sum(numbs, 1);
  arma::vec denbs1 = sum(denbs, 1);
  arma::vec bs = numbs1/denbs1;
  arma::mat Bmat_bs = Bmat.each_col() - bs;
  
  arma::mat num_der_bs_B = 1 / var_bj1T;
  arma::mat der_bs_B = num_der_bs_B;
  der_bs_B.each_col() /= denbs1;
  arma::vec der_ll = arma::zeros(numforms + numvar - 1);
  
  for (s = 1; s < numforms; s++)
  {
    arma::mat der_bs_Br = repmat(der_bs_B.col(s), 1, numforms);
    arma::mat der_llb_jt_Bs =  - 2 * (bj1T_res + Bmat_bs)/var_bj1T % der_bs_Br;
    der_llb_jt_Bs.replace(arma::datum::nan, 0);
    double der_ll_Bs = accu(der_llb_jt_Bs);
    arma::vec tobeadded1 = (bj1T_res.col(s) - bs + B[s]) / var_bj1T.col(s);
    der_ll_Bs += 2 * accu(tobeadded1.elem(arma::find_finite(tobeadded1)));
    der_ll[s - 1] = der_ll_Bs;
  }
  
  arma::mat resb_jt = (bj1T_res + Bmat_bs) / var_bj1T;
  
  if (numvar > 0)
  {
    for (i = 0; i < numvar; i++)
    {
      arma::mat Xmat_i = Xmat(i);
      arma::mat resbX =  - 2 * resb_jt % Xmat_i;
      der_ll[numforms - 1 + i] = accu(resbX.elem(arma::find_finite(resbX)));
    }
  }
  return der_ll;
}

