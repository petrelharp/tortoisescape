all : report-methods.html report-methods.docx


%.html : %.Rmd
	R -e "require(rmarkdown); render(\"$<\")"

%.docx : %.Rmd
	R -e "require(rmarkdown); render(\"$<\",output_format='word_document')"
