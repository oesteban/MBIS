# Makefile to create PDF sub-figures from LaTeX-Files 
# Needed software packages: pdflatex, rubber 
# License: No copyright, just do what the heck you want with it
# vim: set tabstop=4 shiftwidth=4 noexpandtab

fig03-hard_segmentation.pdf: fig03-hard_segmentation.tex
	pdflatex -shell-escape $<
	grep Rerun fig03-hard_segmentation.log && pdflatex $<
	grep Rerun fig03-hard_segmentation.log && pdflatex $<


