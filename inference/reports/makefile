
.PHONY : htmls what

htmls : $(patsubst %.Rmd,%.html,$(wildcard *.Rmd))

%.html : %.Rmd
	R -e 'library(knitr);knit2html("$<")'

model-uncertainty.html : cache/model-description-saved-state.RData

cache/model-description-saved-state.RData : model-description.html

what : 
	echo "$(patsubst %.Rmd,%.html,$(wildcard *.Rmd))"
