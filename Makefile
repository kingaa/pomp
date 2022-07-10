default:
	$(MAKE) -C vignettes default fresh

clean:
	$(MAKE) -C vignettes clean

fresh:
	$(MAKE) -C vignettes fresh
	$(RM) -r _site
