## default templates for pomp's own C snippets.
## This is used in 'pomp.R' and 'builder.R'.

## BASIC TEMPLATES

pomp_templates <- list(
  define=r"{
#define {%variable%} ({%cref%})}",
  undefine=r"{
#undef {%variable%}}",
  file=list(
    header=r"{
/* pomp C snippet file: {%name%} */
/* Time: {%timestamp%} */
/* Salt: {%salt%} */

#include <{%pompheader%}>
#include <R_ext/Rdynload.h>}"
  ),
  utilities=list(
    bspline_basis=list(
      trigger=r"{(?<!periodic_)bspline_basis_eval(?!_deriv)}",
      header=r"{
static bspline_basis_eval_t *__pomp_bspline_basis_eval;
#define bspline_basis_eval(X,K,P,N,Y)  (__pomp_bspline_basis_eval((X),(K),(P),(N),(Y)))}",
      reg=r"{
  __pomp_bspline_basis_eval = (bspline_basis_eval_t *) R_GetCCallable("pomp","bspline_basis_eval");}"
    ),
    bspline_basis_deriv=list(
      trigger=r"{(?<!periodic_)bspline_basis_eval_deriv}",
      header=r"{
static bspline_basis_eval_deriv_t *__pomp_bspline_basis_eval_deriv;
#define bspline_basis_eval_deriv(X,K,P,N,D,Y)  (__pomp_bspline_basis_eval_deriv((X),(K),(P),(N),(D),(Y)))}",
      reg=r"{
  __pomp_bspline_basis_eval_deriv = (bspline_basis_eval_deriv_t *) R_GetCCallable("pomp","bspline_basis_eval_deriv");}"
    ),
    periodic_bspline_basis=list(
      trigger="periodic_bspline_basis_eval(?!_deriv)",
      header=r"{
static periodic_bspline_basis_eval_t *__pomp_periodic_bspline_basis_eval;
#define periodic_bspline_basis_eval(X,T,P,N,Y)  (__pomp_periodic_bspline_basis_eval((X),(T),(P),(N),(Y)))}",
      reg=r"{
  __pomp_periodic_bspline_basis_eval = (periodic_bspline_basis_eval_t *) R_GetCCallable("pomp","periodic_bspline_basis_eval");}"
    ),
    periodic_bspline_basis_deriv=list(
      trigger="periodic_bspline_basis_eval_deriv",
      header=r"{
static periodic_bspline_basis_eval_deriv_t *__pomp_periodic_bspline_basis_eval_deriv;
#define periodic_bspline_basis_eval_deriv(X,T,P,N,D,Y)  (__pomp_periodic_bspline_basis_eval_deriv((X),(T),(P),(N),(D),(Y)))}",
      reg=r"{
  __pomp_periodic_bspline_basis_eval_deriv = (periodic_bspline_basis_eval_deriv_t *) R_GetCCallable("pomp","periodic_bspline_basis_eval_deriv");}"
    ),
    get_userdata_int=list(
      trigger="get_userdata_int",
      header=r"{
static get_userdata_int_t *__pomp_get_userdata_int;
#define get_userdata_int(X)  (__pomp_get_userdata_int(X))}",
      reg=r"{
  __pomp_get_userdata_int = (get_userdata_t *) R_GetCCallable("pomp","get_userdata_int");}"
    ),
    get_userdata_double=list(
      trigger="get_userdata_double",
      header=r"{
static get_userdata_double_t *__pomp_get_userdata_double;
#define get_userdata_double(X)  (__pomp_get_userdata_double(X))}",
      reg=r"{
  __pomp_get_userdata_double = (get_userdata_double_t *) R_GetCCallable("pomp","get_userdata_double");}"
    ),
    get_userdata=list(
      trigger=r"{get_userdata(?!_)}",
      header=r"{
static get_userdata_t *__pomp_get_userdata;
#define get_userdata(X)  (__pomp_get_userdata(X))}",
      reg=r"{
  __pomp_get_userdata = (get_userdata_t *) R_GetCCallable("pomp","get_userdata");}"
    )
  ),
  stackhandling=r"{
static int __pomp_load_stack = 0;
void __pomp_load_stack_incr (void) {++__pomp_load_stack;}
void __pomp_load_stack_decr (int *val) {*val = --__pomp_load_stack;}}",
  registration=list(
    header=r"{
void R_init_{%name%} (DllInfo *info) {}",
    main=r"(
  R_RegisterCCallable("{%name%}", "{%fun%}", (DL_FUNC) {%fun%});)",
    footer="\n}\n"
  )
)

## WORKHORSE-SPECIFIC TEMPLATES

workhorse_templates <- list(
  rinit=list(
    slotname="rinit",
    Cname="__pomp_rinit",
    proto=quote(rinit(...)),
    header=r"{
void __pomp_rinit (double *__x, const double *__p, double t,
   const int *__stateindex, const int *__parindex, const int *__covindex,
   const double *__covars) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      )
    )
  ),
  dinit=list(
    slotname="dinit",
    Cname="__pomp_dinit",
    proto=quote(dinit(...)),
    header=r"{
void __pomp_dinit (double *__loglik, const double *__x, const double *__p,
  double t, const int *__stateindex, const int *__parindex, const int *__covindex,
  const double *__covars) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      loglik=list(
        names=r"{loglik}",
        cref=r"{__loglik[0]}"
      )
    )
  ),
  rmeasure=list(
    slotname="rmeasure",
    Cname="__pomp_rmeasure",
    proto=quote(rmeasure(...)),
    header=r"{
void __pomp_rmeasure (double *__y, const double *__x, const double *__p,
  const int *__obsindex, const int *__stateindex, const int *__parindex,
  const int *__covindex, const double *__covars, double t) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      obs=list(
        names=quote(obsnames),
        cref=r"{__y[__obsindex[{%v%}]]}"
      )
    )
  ),
  emeasure=list(
    slotname="emeasure",
    Cname="__pomp_emeasure",
    proto=quote(emeasure(...)),
    header=r"{
void __pomp_emeasure (double *__f, const double *__x, const double *__p,
  const int *__obsindex, const int *__stateindex, const int *__parindex,
  const int *__covindex, const double *__covars, double t) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      expectations=list(
        names=quote(paste0("E_",obsnames)),
        cref=r"{__f[__obsindex[{%v%}]]}"
      )
    )
  ),
  vmeasure=list(
    slotname="vmeasure",
    Cname="__pomp_vmeasure",
    proto=quote(vmeasure(...)),
    header=r"{
void __pomp_vmeasure (double *__f, const double *__x, const double *__p,
  const int *__vmatindex, const int *__stateindex, const int *__parindex,
  const int *__covindex, const double *__covars, int nobs, double t) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      variance_matrix=list(
        names=quote(paste0("V_",outer(obsnames,obsnames,paste,sep="_"))),
        cref=r"{__f[__vmatindex[{%v%}]]}"
      )
    )
  ),
  dmeasure=list(
    slotname="dmeasure",
    Cname= "__pomp_dmeasure",
    proto=quote(dmeasure(log,...)),
    header=r"{
void __pomp_dmeasure (double *__lik, const double *__y, const double *__x,
  const double *__p, int give_log, const int *__obsindex,
  const int *__stateindex, const int *__parindex, const int *__covindex,
  const double *__covars, double t) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      obs=list(
        names=quote(obsnames),
        cref=r"{__y[__obsindex[{%v%}]]}"
      ),
      lik=list(
        names="lik",
        cref=r"{__lik[0]}"
      )
    )
  ),
  step.fn=list(
    slotname="step.fun",
    Cname="__pomp_stepfn",
    proto=quote(step.fun(...)),
    header=r"{
void __pomp_stepfn (double *__x, const double *__p, const int *__stateindex,
  const int *__parindex, const int *__covindex, const double *__covars,
  double t, double dt) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      )
    )
  ),
  rate.fn=list(
    slotname="rate.fun",
    Cname="__pomp_ratefn",
    proto=quote(rate.fun(j,...)),
    header=r"{
double __pomp_ratefn (int j, double t, double *__x, const double *__p,
  const int *__stateindex, const int *__parindex, const int *__covindex,
  const double *__covars) {
  double rate = 0.0;}",
    footer=r"{
  return rate;
}}",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      )
    )
  ),
  dprocess=list(
    slotname="dprocess",
    Cname="__pomp_dproc",
    proto=quote(dprocess(...)),
    header=r"{
void __pomp_dproc (double *__loglik, const double *__x1, const double *__x2,
  double t_1, double t_2, const double *__p, const int *__stateindex,
  const int *__parindex, const int *__covindex, const double *__covars) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      before=list(
        names=quote(paste0(statenames,"_1")),
        cref=r"{__x1[__stateindex[{%v%}]]}"
      ),
      after=list(
        names=quote(paste0(statenames,"_2")),
        cref=r"{__x2[__stateindex[{%v%}]]}"
      ),
      loglik=list(
        names="loglik",
        cref=r"{__loglik[0]}"
      )
    )
  ),
  skeleton=list(
    slotname="skeleton",
    Cname="__pomp_skelfn",
    proto=quote(skeleton(...)),
    header=r"{
void __pomp_skelfn (double *__f, const double *__x, const double *__p,
  const int *__stateindex, const int *__parindex, const int *__covindex,
  const double *__covars, double t) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      states=list(
        names=quote(statenames),
        cref=r"{__x[__stateindex[{%v%}]]}"
      ),
      derivs=list(
        names=quote(paste0("D",statenames)),
        cref=r"{__f[__stateindex[{%v%}]]}"
      )
    )
  ),
  fromEst=list(
    slotname="fromEst",
    Cname="__pomp_from_trans",
    proto=quote(from.trans(...)),
    header=r"{
void __pomp_from_trans (double *__p, const double *__pt, const int *__parindex) {}",
    footer="\n}\n",
    vars=list(
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      transforms=list(
        names=quote(paste0("T_",paramnames)),
        cref=r"{__pt[__parindex[{%v%}]]}"
      )
    )
  ),
  toEst=list(
    slotname="toEst",
    Cname="__pomp_to_trans",
    proto=quote(to.trans(...)),
    header=r"{
void __pomp_to_trans (double *__pt, const double *__p, const int *__parindex) {}",
    footer="\n}\n",
    vars=list(
      covars=list(
        names=quote(covarnames),
        cref=r"{__covars[__covindex[{%v%}]]}"
      ),
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      transforms=list(
        names=quote(paste0("T_",paramnames)),
        cref=r"{__pt[__parindex[{%v%}]]}"
      )
    )
  ),
  rprior=list(
    slotname="rprior",
    Cname="__pomp_rprior",
    proto=quote(rprior(...)),
    header=r"{
void __pomp_rprior (double *__p, const int *__parindex) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      )
    )
  ),
  dprior=list(
    slotname="dprior",
    Cname="__pomp_dprior",
    proto=quote(dprior(log,...)),
    header=r"{
void __pomp_dprior (double *__lik, const double *__p, int give_log,
  const int *__parindex) {}",
    footer="\n}\n",
    vars=list(
      params=list(
        names=quote(paramnames),
        cref=r"{__p[__parindex[{%v%}]]}"
      ),
      lik=list(
        names="lik",
        cref=r"{__lik[0]}"
      )
    )
  )
)
