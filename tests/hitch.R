options(digits=3)

library(pomp2)

try(hitch())

hitch(templates=pomp2:::workhorse_templates)

try(hitch(templates=pomp2:::workhorse_templates,
  statenames=c("a","b"),paramnames=c("b","c")))

hitch(step.fn=Csnippet("int bob; bob = 3"),
  templates=pomp2:::workhorse_templates,
  compile=FALSE,cfile="bob") %>% names()

