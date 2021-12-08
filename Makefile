texwfiles=$(wildcard *.texw)
pyfiles=$(patsubst %.texw,%.py,$(texwfiles))
pdffiles=$(patsubst %.texw,%.pdf,$(texwfiles))

all: $(pyfiles)

pdfs: $(pdffiles)


%.py: %.texw
	ptangle $<

%.tex: %.texw
	pweave -f texpygments $<

%.pdf: %.tex pygments.sty
	pdflatex $<

pygments.sty: Makefile
	pygmentize -f tex -S friendly > pygments.sty
	# pygmentize -f tex -S tango > pygments.sty
	# pygmentize -f tex -S lovelace > pygments.sty

