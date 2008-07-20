## draw a set of Np particles from the user-specified distribution
setMethod(
          "particles",
          "mif",
          function (object, Np = 1, center = coef(object), sd = 0, ...) {
            if ((length(sd)==1) && (sd == 0)) {
              sd <- rep(0,length(center))
              names(sd) <- names(center)
            }
            if (is.null(names(center)) || is.null(names(sd)))
              stop("particles error: 'center' and 'sd' must have names",call.=FALSE)
            if (length(sd)!=length(center))
              stop("particles error: 'center' and 'sd' must be of equal length",call.=FALSE)
            x <- try(
                     do.call(
                             object@particles,
                             c(
                               list(Np=Np,center=center,sd=sd),
                               object@userdata
                               )
                             ),
                     silent=FALSE
                     )
            if (inherits(x,'try-error'))
              stop("particles error: error in user-specified 'particles' function",call.=FALSE)
            if (
                !is.matrix(x) ||
                Np!=ncol(x) ||
                is.null(rownames(x))
                )
              stop("particles error: user 'particles' function must return a matrix with Np columns and rownames",call.=FALSE)
            x
          }
          )
