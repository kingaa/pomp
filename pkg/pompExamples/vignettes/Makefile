RSCRIPT = Rscript --vanilla
RM = rm -f

vignettes: bsflu.html
	$(RM) -r figure

%.html: %.Rmd
	$(RSCRIPT) -e 'knitr::knit("$*.Rmd")'
	$(RSCRIPT) -e 'markdown::markdownToHTML("$*.md","$*.html")'
	$(RM) $*.md
