texwfiles=$(wildcard *.texw)
pyfiles=$(patsubst %.texw,%.py,$(texwfiles))
pdffiles=$(patsubst %.texw,%.pdf,$(texwfiles))

all: pdfs pys

pys: $(pyfiles)

pdfs: $(pdffiles)

%.py: %.texw
	ptangle $<

%.tex: %.texw
	pweave -f texminted $<

%.pdf: %.tex
	pdflatex -shell-escape $<
	pdflatex -shell-escape $<

