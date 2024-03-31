TARBALL = $(wildcard *.tar.gz)

rchk: rchk.sif
	singularity run rchk.sif $(PWD)/$(TARBALL) | tee rchk.out

rchk.sif:
	singularity pull rchk.sif shub://kalibera/rchk:def

clean:
	$(RM) $(TARBALL)
	$(RM) -r build lib libsonly

fresh: clean
	$(RM) rchk.sif
