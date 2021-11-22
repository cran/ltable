
#include <Rcpp.h>
#include <RcppGSL.h> 
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
using namespace Rcpp;

struct data {
  size_t n;
  size_t p;
  gsl_matrix * xdata;
  gsl_vector * y;
};

int floglin(const gsl_vector * beta, void *data,
            gsl_vector * f);
int dfloglin (const gsl_vector * beta, void *data,
              gsl_matrix * J);
void callback(const size_t iter, void *params,
              const gsl_multifit_nlinear_workspace *w);


int floglin(const gsl_vector * beta, void *data,
            gsl_vector * f) {
  
  size_t n = ((struct data *)data)->n;
  size_t p = ((struct data *)data)->p;
  gsl_matrix * xdata = ((struct data *)data)->xdata;
  gsl_vector * y = ((struct data *)data)->y;
  
  gsl_matrix * LP = gsl_matrix_alloc (n, p);
  gsl_matrix_memcpy(LP, xdata) ;
 
  size_t i, j;
  for (i = 0; i < n; i++)
  {
    gsl_vector_view z = gsl_matrix_row(LP,i);
    gsl_vector_mul( &z.vector, beta);
    double sum_v = 0;
    for(j=0; j<p; j++) {sum_v += gsl_vector_get(&z.vector,j);}
    double Yi = exp(sum_v);

    gsl_vector_set (f, i, Yi - gsl_vector_get(y,i));
  }
  gsl_matrix_free(LP);
  
  return GSL_SUCCESS;
}
  int dfloglin (const gsl_vector * beta, void *data,
             gsl_matrix * J) {
    
    size_t n = ((struct data *)data)->n;
    size_t p = ((struct data *)data)->p;
    gsl_matrix * xdata = ((struct data *)data)->xdata;
    gsl_matrix * LP = gsl_matrix_alloc (n, p);
    gsl_matrix_memcpy(LP, xdata) ;

    size_t i;
    size_t j;
  
    for (i = 0; i < n; i++)
    {
      gsl_vector_view z = gsl_matrix_row(LP,i);
      gsl_vector_mul( &z.vector, beta);
      double sum_v = 0;
      for(j=0; j<p; j++) {sum_v += gsl_vector_get(&z.vector,j);}
      double Basei=exp(sum_v);
      for (j = 0; j < p; j++)
      /* Jacobian matrix J(i,j) = dfi / dxj, */
      gsl_matrix_set (J, i, j, gsl_matrix_get(xdata,i,j)*Basei);
    }
    gsl_matrix_free(LP);
    return GSL_SUCCESS;
  }
    
  
  
  void callback(const size_t iter, void *params,
           const gsl_multifit_nlinear_workspace *w)
  {
    /*    gsl_vector *f = gsl_multifit_nlinear_residual(w);
    gsl_vector *beta = gsl_multifit_nlinear_position(w);
    double rcond;
    compute reciprocal condition number of J(x) 
    gsl_multifit_nlinear_rcond(&rcond, w);*/
  }

//[[Rcpp::export]] 



Rcpp::List gsl_poisson(const RcppGSL::matrix<double> & G) {
 
 
   const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
   gsl_multifit_nlinear_workspace *w;
   gsl_multifit_nlinear_fdf fdf;
   gsl_multifit_nlinear_parameters fdf_params =
   gsl_multifit_nlinear_default_parameters();
   
   const size_t n = G -> size1;
   const size_t k = G -> size2;
   const size_t p = k-1;
   gsl_matrix * DM = gsl_matrix_alloc (n, k);
   gsl_matrix_memcpy(DM, G) ;
   
   gsl_vector_view yview= gsl_matrix_column(DM, 0);
   gsl_matrix_view Xview = gsl_matrix_submatrix(DM,0,1,n,p);
   gsl_matrix * xdata=&Xview.matrix;
   gsl_vector * y=&yview.vector;
   
  struct data d = { n, p, xdata, y };  
   
  gsl_vector *f;
  gsl_matrix *J;
  gsl_matrix *covar = gsl_matrix_alloc (p, p);
  gsl_vector * beta_init = gsl_vector_alloc (p);
  gsl_vector_set_all(beta_init,1); /* starting values */
  gsl_vector_view beta = gsl_vector_subvector(beta_init, 0, p);
  
  /* WEIGHTS */
  gsl_vector * rowdata = gsl_vector_alloc (p); 
  gsl_vector * weights = gsl_vector_alloc (n); 
  size_t i, j;
  
  for (i = 0; i < n; i++)
  {
    gsl_vector_view z = gsl_matrix_row(xdata,i);
    gsl_vector_memcpy(rowdata, &z.vector) ;
    gsl_vector_mul( rowdata, &beta.vector);
    double sum_v = 0;
    for(j=0; j<p; j++) {sum_v += gsl_vector_get(&z.vector,j);}  
    double mui=exp(sum_v);   
  gsl_vector_set (weights, i, 1.0 / abs(mui-gsl_vector_get(y,i)));
  };
 
 // gsl_rng * r;
  double chisq, chisq0;
  int status, info;
  
  
  const double xtol = 1e-8;
  const double gtol = 1e-8;
  const double ftol = 0.0;
//  gsl_rng_env_setup();
//  r = gsl_rng_alloc(gsl_rng_default);

/* define the function to be minimized */
  fdf.f = floglin;
  fdf.df = dfloglin;
  fdf.fvv = NULL;
  fdf.n = n;
  fdf.p = p;
  fdf.params = &d;

  /* allocate workspace with default parameters */
  
  w = gsl_multifit_nlinear_alloc (T, &fdf_params, n, p);
   

  /* initialize solver with starting point and weights */ 
  gsl_multifit_nlinear_winit(&beta.vector, weights, &fdf, w);
  
  /* compute initial cost function */
 f = gsl_multifit_nlinear_residual(w);
 gsl_blas_ddot(f, f, &chisq0);
  
  /* solve the system with a maximum of 100 iterations */
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol,
                                       callback, NULL, &info, w);
  
   /* compute covariance of best fit parameters */
  J = gsl_multifit_nlinear_jac(w);
  gsl_multifit_nlinear_covar (J, 0.0, covar);
  
   /* compute final cost */
  gsl_blas_ddot(f, f, &chisq);
 
  
  gsl_vector_view betavars=gsl_matrix_diagonal(covar);
  
  Rcpp:: NumericVector betas(p), errors (p);
  
  for (i = 0; i < p; i++)
  {
    betas[i]=gsl_vector_get(w->x, i);
    errors[i]=sqrt(gsl_vector_get(&betavars.vector, i));
  }
  double rcond;
  /* compute reciprocal condition number of J(x) */
  gsl_multifit_nlinear_rcond(&rcond, w);
 
 
  Rcpp::List res = Rcpp::List::create( _["betas"] = betas, 
                                       _["errors"]=errors,
                                       _["initial |f(x)|"]=sqrt(chisq0),
                                       _["final |f(x)|"]=sqrt(chisq),
                                       _["Jacobian reciprocal condition number"]=1.0 / rcond,
                                       _["number of iterations"]=gsl_multifit_nlinear_niter(w),
                                       _["reason for stopping"]=(info == 1) ? "small step size" : "small gradient",
                                       _["chisq/dof"]= chisq / (n-p),
                                       _["status"]=gsl_strerror (status) );
  gsl_multifit_nlinear_free (w);
  gsl_matrix_free (covar);
  gsl_matrix_free(DM);
  gsl_vector_free(rowdata);
  gsl_vector_free(weights);
  
  

 return res;  
  }
  
