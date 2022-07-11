default:
	$(MAKE) -C vignettes default fresh

serve:
	jekyll serve

clean:
	$(MAKE) -C vignettes clean

fresh:
	$(MAKE) -C vignettes fresh
	$(RM) -r _site
