png(filename="mcap-%02d.png",res=100)

library(pomp)
library(dplyr)
set.seed(722855899)

profile_design(
  theta=seq(-2,2,length=100),
  lower=c(a=0,b=0),
  upper=c(a=1,b=2),
  nprof=10
) %>%
  mutate(
    logLik=rnorm(
      n=length(theta),
      mean=-theta^2-0.1*(theta+0.5)^4-300,
      sd=0.3*(a+b)
    )
  ) -> x

mcp <- mcap(x$logLik,x$theta)

plot(logLik~theta,data=x,pch=16)
lines(smoothed~parameter,data=mcp$fit,col=4,lwd=3)
lines(quadratic~parameter,data=mcp$fit,col=7,lwd=3)
abline(v=mcp$ci,lwd=3,lty=2,col=1)
abline(v=mcp$mle,lwd=3,lty=3,col=1)
text(x=-0.5,y=-305,labels=bquote(MLE==.(signif(mcp$mle,2))))
text(x=-0.5,y=-306,labels=bquote(se[stat]==.(signif(mcp$se_stat,3))))
text(x=-0.5,y=-306.5,labels=bquote(se[mc]==.(signif(mcp$se_mc,3))))
text(x=-0.5,y=-307,labels=bquote(se[total]==.(signif(mcp$se,3))))

dev.off()
