# for recording how stuff got made:
#


tort_272_info/pcs.csv : tort_272_info/covariance.csv
	Rscript ../covmat/write-PCs.R covariance.csv pcs.csv

tort_276_info/covariance.csv : tort_276_info/276torts.covar
	# did some stuff in R here
	
tort_276_info/all_angsd_snps.pwp.csv : pairwisePi/allSNPs.pwp
	# did some stuff in R here

tort_276_info/all_sites.pwp.csv : pairwisePi/allllllSites.pwp
	# did some stuff in R here
