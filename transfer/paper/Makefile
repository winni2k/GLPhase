### Makefile --- converge paper

## Author: wkretzsch@gmail.com
## Version: $Id: Makefile,v 0.0 2013/05/20 16:49:20 winni Exp $
## Keywords: 
## X-URL: 

output = .output
backup = .backup
base   = transfer

all:
	-mkdir $(output)
	-cp *.sty *.tex *.bib $(output)
	-mv *.pdf *.eps $(output)
	cd $(output) && ln -s ../img img; latexmk -pdf -pv $(base).tex && cp $(base).pdf ../

clean:
	-rm $(output)/*

backup:
	-mkdir $(backup)
	cp  $(output)/$(base).pdf $(backup)/$(base).pdf

### Makefile ends here
