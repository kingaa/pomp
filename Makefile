default:
	$(MAKE) -C vignettes default fresh

serve:
	jekyll serve --drafts -b /pomp

clean:
	$(MAKE) -C vignettes clean

fresh:
	$(MAKE) -C vignettes fresh
	$(RM) -r _site .jekyll-cache
