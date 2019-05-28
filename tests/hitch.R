options(digits=3)

library(pomp)

try(hitch())

hitch(templates=pomp:::workhorse_templates)

try(hitch(templates=pomp:::workhorse_templates,
  statenames=c("a","b"),paramnames=c("b","c")))

hitch(step.fn=Csnippet("int bob; bob = 3"),
  templates=pomp:::workhorse_templates,
  compile=FALSE,cfile="bob") %>% names()

