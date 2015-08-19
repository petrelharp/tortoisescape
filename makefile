# for recording how stuff got made:
#

## PCs, all sites

tort_272_info/allSites_pcs-PC-maps.pdf tort_272_info/allSites_pcs-maps-with-PCs.pdf : tort_272_info/allSites_pcs.csv
	cd tort_272_info; Rscript  ../visualization/plot-pcs.R allSites_pcs.csv

tort_272_info/allSites_pcs.csv : tort_272_info/allSites_covar.csv tort_272_info/weights-nodups.csv
	Rscript visualization/write-PCs.R $< $@ tort_272_info/weights-nodups.csv

## PCs, polymorphic sites

tort_272_info/pcs-PC-maps.pdf tort_272_info/pcs-maps-with-PCs.pdf : tort_272_info/pcs.csv
	cd tort_272_info; Rscript  ../visualization/plot-pcs.R pcs.csv

tort_272_info/pcs.csv : tort_272_info/covariance.csv tort_272_info/weights-nodups.csv
	Rscript visualization/write-PCs.R tort_272_info/covariance.csv $@ tort_272_info/weights-nodups.csv

tort_272_info/weighted-pcs.csv : tort_272_info/covariance.csv tort_272_info/weights.csv
	Rscript visualization/write-PCs.R tort_272_info/covariance.csv $@ tort_272_info/weights.csv

tort_272_info/weights.csv tort_272_info/weights-nodups.csv : 
	cd tort_272_info; Rscript make-weights.R

tort_272_info/covariance.csv : covmat/allSNPs-cpp_labeled.txt
	# did some stuff in R here 2/16/2015:
	#   added etort IDs

tort_272_info/all_angsd_snps.pwp.csv : pairwisePi/allSNPs.pwp
	# did some stuff in R here

tort_272_info/all_sites.pwp.csv : pairwisePi/allllllSites.pwp
	# did some stuff in R here

tort_272_info/all_sites_minus_het.pwp.csv : pairwisePi/allllllSites.pwp
	# subtracted average of heterozygosities
	# and used df.to.sym and sym.to.df in inference/input-output-fns.R

tort_272_info/all_sites_pwp_pngs : tort_272_info/geog_distance.csv tort_272_info/all_sites.pwp.csv
	cd tort_272_info; Rscript ../visualization/plot-pwp-pngs.R geog_distance.csv all_sites.pwp.csv . all_sites_pwp_pngs

tort_272_info/all_angsd_snps_pngs : tort_272_info/geog_distance.csv tort_272_info/all_angsd_snps.pwp.csv
	cd tort_272_info; Rscript ../visualization/plot-pwp-pngs.R geog_distance.csv all_angsd_snps.pwp.csv . all_angsd_snps_pngs

tort_272_info/all_sites_pwp_minus_het_pngs : tort_272_info/geog_distance.csv tort_272_info/all_sites_minus_het.pwp.csv
	cd tort_272_info; Rscript ../visualization/plot-pwp-pngs.R geog_distance.csv all_sites_minus_het.pwp.csv . all_sites_pwp_minus_het_pngs

tort_272_info/all_angsd_snps_distance_maps : tort_272_info/all_angsd_snps.pwp.csv
	cd tort_272_info; Rscript ../visualization/distance-maps.R geog_distance.csv all_angsd_snps.pwp.csv . all_angsd_snps_distance_maps

%_distance_maps : %.pwp.csv
	cd tort_272_info; Rscript ../visualization/distance-maps.R geog_distance.csv $(notdir $<) . $(notdir $@)
	make $@/index.html

%/index.html :
	R --vanilla -e "dirname=\"../$*\";knitr::knit(\"visualization/index.Rmd\",output=\"$*/index.md\")"
	pandoc -s $*/index.md -t slidy -o $@
