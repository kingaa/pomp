## default templates for pomp's own C snippets.
## This is used in 'pomp.R', 'plugins.R', and 'builder.R'.
snippet_templates <- list(
  initializer=list(
    Cname="__pomp_rinit",
    needs=c("paramnames","statenames","covarnames"),
    header="\nvoid __pomp_rinit (double *__x, const double *__p, double t, const int *__stateindex, const int *__parindex, const int *__covindex, const double *__covars)\n{\n",
    footer="\n}\n\n"
  ),
  rmeasure=list(
    Cname="__pomp_rmeasure",
    needs=c("paramnames","statenames","obsnames","covarnames"),
    header="\nvoid __pomp_rmeasure (double *__y, const double *__x, const double *__p, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  dmeasure=list(
    Cname= "__pomp_dmeasure",
    needs=c("paramnames","statenames","obsnames","covarnames","lik"),
    header="\nvoid __pomp_dmeasure (double *__lik, const double *__y, const double *__x, const double *__p, int give_log, const int *__obsindex, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  step.fn=list(
    Cname="__pomp_stepfn",
    needs=c("paramnames","covarnames","statenames"),
    header="\nvoid __pomp_stepfn (double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars, double t, double dt)\n{\n",
    footer="\n}\n\n"
  ),
  rate.fn=list(
    Cname="__pomp_ratefn",
    needs=c("paramnames","covarnames","statenames"),
    header="\ndouble __pomp_ratefn (int j, double t, double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __covdim, const double *__covars)\n{\n  double rate = 0.0;  \n",
    footer="  return rate;\n}\n\n"
  ),
  dens.fn=list(
    Cname="__pomp_densfn",
    needs=c("paramnames","covarnames","before_n_after","loglik"),
    header="\nvoid __pomp_densfn (double *__loglik, const double *__x1, const double *__x2, double t_1, double t_2, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars)\n{\n",
    footer="\n}\n\n"
  ),
  skeleton=list(
    Cname="__pomp_skelfn",
    needs=c("paramnames","covarnames","statenames","derivs"),
    header="\nvoid __pomp_skelfn (double *__f, const double *__x, const double *__p, const int *__stateindex, const int *__parindex, const int *__covindex, int __ncovars, const double *__covars, double t)\n{\n",
    footer="\n}\n\n"
  ),
  fromEstimationScale=list(
    Cname="__pomp_from_trans",
    needs=c("paramnames","transforms","covarnames"),
    header="\nvoid __pomp_from_trans (double *__pt, const double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  toEstimationScale=list(
    Cname="__pomp_to_untrans",
    needs=c("paramnames","transforms","covarnames"),
    header="\nvoid __pomp_to_untrans (double *__pt, const double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  rprior=list(
    Cname="__pomp_rprior",
    needs=c("paramnames","covarnames"),
    header="\nvoid __pomp_rprior (double *__p, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  ),
  dprior=list(
    Cname="__pomp_dprior",
    needs=c("paramnames","covarnames","lik"),
    header="\nvoid __pomp_dprior (double *__lik, const double *__p, int give_log, const int *__parindex)\n{\n",
    footer="\n}\n\n"
  )
)
