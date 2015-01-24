# for recording how stuff got made:
#


tort_272_info/pcs.csv : tort_272_info/covariance.csv
	Rscript covmat/write-PCs.R tort_272_info/covariance.csv tort_272_info/pcs.csv

tort_272_info/covariance.csv : tort_276_info/276torts.covar
	# did some stuff in R here

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
