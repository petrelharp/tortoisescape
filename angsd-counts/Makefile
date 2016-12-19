%.perindividual.hists : %.bin
	Rscript ../count-utils/get-coverage-histogram.R $< 100 > $@

%.persite.stats : %.bin
	Rscript ../count-utils/get-coverage-stats.R $< > $@
