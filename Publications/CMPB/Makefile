# Makefile to create PDF documents from LaTeX-Files 
# Needed software packages: pdflatex, rubber 
# License: No copyright, just do what the heck you want with it
# vim: set tabstop=4 shiftwidth=4 noexpandtab
# Partially from: http://blog.dbrgn.ch/2011/2/22/a-simple-latex-makefile/

SUBDIRS = figures 

all: COMM_3778.pdf

package: COMM_3778.zip

COMM_3778.zip: COMM_3778/COMM_3778.pdf
	rm -f COMM_3778/COMM_3778.{aux,log,pdf,blg};\
	zip -r COMM_3778.zip COMM_3778/

COMM_3778/COMM_3778.pdf: COMM_3778/COMM_3778.tex COMM_3778/COMM_3778.bbl COMM_3778/COMM_3778.aux
	cd COMM_3778/; \
	pdflatex -shell-escape COMM_3778.tex; \
	grep Rerun COMM_3778.log && pdflatex COMM_3778.tex; \
	grep Rerun COMM_3778.log && pdflatex COMM_3778.tex ; \
	touch COMM_3778.bbl ;
	touch $@ ;


COMM_3778/COMM_3778.bbl: COMM_3778/COMM_3778.aux COMM_3778/COMM_3778.tex
	cd COMM_3778; \
	bibtex COMM_3778

COMM_3778/COMM_3778.aux: COMM_3778/COMM_3778.tex expanded-figures
	cd COMM_3778; \
	pdflatex -shell-escape COMM_3778.tex

expanded-figures:
	for fname in $$( sed -n "s/^.*\\includegraphics\[.*\]{\(.*\)}/\\1/p" COMM_3778/COMM_3778.tex ); do \
		convert figures/$${fname}.pdf -density 300 COMM_3778/$${fname}.pdf; \
	done

COMM_3778/COMM_3778.tex: *.tex COMM_3778.bib
	mkdir -p COMM_3778
	cp COMM_3778.bib COMM_3778/
	cp COMM_3778.bst COMM_3778/
	latexpand COMM_3778.tex | sed '/\graphicspath/d' > COMM_3778/COMM_3778.tex

COMM_3778.pdf: COMM_3778.tex COMM_3778.bbl COMM_3778.aux $(SUBDIRS)
	pdflatex -shell-escape $< 
	grep Rerun COMM_3778.log && pdflatex $< 
	grep Rerun COMM_3778.log && pdflatex $< 
	touch COMM_3778.bbl
	touch $@

COMM_3778.bbl: COMM_3778.aux COMM_3778.bib
	bibtex COMM_3778

COMM_3778.aux: *.tex COMM_3778.bib $(SUBDIRS)
	pdflatex -shell-escape COMM_3778.tex

COMM_3778.bib: MBIS.bib
	cat MBIS.bib > COMM_3778.bib

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

cleanall: clean
	rm -f *.pdf *.acn *.glo *.ist *.toc *.spl *.dvi *.ps *.ttt *.fff *.lot COMM_3778.txt COMM_3778.odp COMM_3778.bib COMM_3778.zip
	rm -rf COMM_3778/

clean:
	rm -f *.aux *.log *.out *.blg *.bbl package


.PHONY: all pdf clean cleanall 
.SILENT: pdf

