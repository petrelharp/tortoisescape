paste <(ls run_*/*/model.score) <(cat run_*/*/model.score) | sort -k 2 -n > all.model.scores
