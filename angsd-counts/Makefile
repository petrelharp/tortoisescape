%.perindividual.hists : %.bin
	Rscript ../../count-utils/get-coverage-histogram.R $< 30 > $@

%.persite.stats : %.bin
	Rscript ../../count-utils/get-coverage-stats.R $< > $@
