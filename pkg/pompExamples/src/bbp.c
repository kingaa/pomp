/* pomp model file: _bombay_plague */

#include <pomp.h>
#include <R_ext/Rdynload.h>

 
#define Beta	(__p[__parindex[0]])
#define delta	(__p[__parindex[1]])
#define mu	(__p[__parindex[2]])
#define gamma	(__p[__parindex[3]])
#define sigma	(__p[__parindex[4]])
#define theta	(__p[__parindex[5]])
#define ratio	(__p[__parindex[6]])
#define x	(__x[__stateindex[0]])
#define y	(__x[__stateindex[1]])
#define n	(__x[__stateindex[2]])
#define deaths	(__y[__obsindex[0]])
#define Dx	(__f[__stateindex[0]])
#define Dy	(__f[__stateindex[1]])
#define Dn	(__f[__stateindex[2]])
#define TBeta	(__pt[__parindex[0]])
#define Tdelta	(__pt[__parindex[1]])
#define Tmu	(__pt[__parindex[2]])
#define Tgamma	(__pt[__parindex[3]])
#define Tsigma	(__pt[__parindex[4]])
#define Ttheta	(__pt[__parindex[5]])
#define Tratio	(__pt[__parindex[6]])
#define lik	(__lik[0])

void _bombay_plague_rmeasure (double *__y, double *__x, double *__p, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 deaths=rnbinom_mu(theta,ratio*exp(y)); 
}


void _bombay_plague_dmeasure (double *__lik, double *__y, double *__x, double *__p, int give_log, int *__obsindex, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{
 lik=dnbinom_mu(deaths,theta,ratio*exp(y),give_log); 
}


void _bombay_plague_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)
{

                         double X = exp(x);
                         double Y = exp(y);
                         double dx, dy, dn, dW, ito;
                         dx = (mu*(1.0/X-1)+(delta-Beta)*Y)*dt;
                         dy = (Beta*X+delta*(Y-1)-gamma-mu)*dt;
                         dn = -delta*Y*dt;
                         dW = rnorm(0,sigma*sqrt(dt));
                         ito = 0.5*sigma*sigma*dt;
                         x += dx - Beta*Y*(dW-Beta*Y*ito);
                         y += dy + Beta*X*(dW+Beta*X*ito);
                         n += dn;
                          
}


void _bombay_plague_skelfn (double *__f, double *__x, double *__p, int *__stateindex, int *__parindex, int *__covindex, int __ncovars, double *__covars, double t)
{

                        double X = exp(x);
                        double Y = exp(y);
                        Dx = mu*(1.0/X-1)+(delta-Beta)*Y;
                        Dy = Beta*X+delta*(Y-1)-gamma-mu;
                        Dn = -delta*Y;
                         
}


void _bombay_plague_rprior (double *__p, int *__parindex)
{
   error("'rprior' not defined"); 
}


void _bombay_plague_dprior (double *__lik, double *__p, int give_log, int *__parindex)
{
   error("'dprior' not defined"); 
}

#undef Beta
#undef delta
#undef mu
#undef gamma
#undef sigma
#undef theta
#undef ratio
#undef x
#undef y
#undef n
#undef deaths
#undef Dx
#undef Dy
#undef Dn
#undef TBeta
#undef Tdelta
#undef Tmu
#undef Tgamma
#undef Tsigma
#undef Ttheta
#undef Tratio
