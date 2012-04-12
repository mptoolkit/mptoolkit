# $Id$
#
# makefile for LaTeX.
#
# This extracts figure dependencies (currently recognized only as .eps files, this could be improved).
# Automatically generates .eps files from xmgrace .agr files.
# There are some hacks for fixing the bounding box which xmgrace seems to get wrong.  This probably
# only works in special cases, the bounding box probably needs adjusting (or better still, I think
# there are tools to fix the bounding box automatically?)
# Also, in case of error the .dvi file is adjusted so that the modification time is the same as the
# modification time of the .tex source.  This is done because latex updates the .dvi file even if there
# is an error.  Consider the following scenario:
# 1. The latex compiler encounters an error and shows a prompt.
# 2. The user fixes the error in some editor window.
# 3. The user then stops the latex compiler.  At this point, the modification time of the .dvi file
# is LATER than the modification time of the .tex file, and a subsequent make will fail.
#
# Additions:
# If the .tex file has dependencies other than .eps files, these are not recognized.  It would be good
# to generalize the .figures file to generate the full dependency tree.
# The .pdf target needs reviewing.  It works, but I'm not sure if the semantics of pdflatex 
# are the same as ordinary latex.
# Some more targets for generating images would be nice.  xfig?
# Bibtex integration?  Does anyone use bibtex?
# There are lots of latex makefiles on the internet where ideas for solving these problems can be obtained.
# Example: http://www.acoustics.hut.fi/u/mairas/UltimateLatexMakefile/

# bounding box to use for postscript from grace (no quotes!)
GraceBoundingBox = 50 15 575 750
#GraceBoundingBox = 0 15 575 750

#default : all

# all .tex sources in the current directory
TEXSRC = $(shell egrep -l '^[^%]*\\begin\{document\}' *.tex)

DVIFILES = $(TEXSRC:%.tex=%.dvi)
PSFILES = $(TEXSRC:%.tex=%.ps)
PDFFILES = $(TEXSRC:%.tex=%.pdf)

# default target is to make all .dvi files
all : $(TEXSRC:.tex=.dvi)
.PHONY : all clean ps pdf

# include the figure dependencies
include $(TEXSRC:.tex=.figures)

ps : $(PSFILES)

pdf : $(PDFFILES)

# 'make figures' generates .eps for all figures
figures : $(TEXSRC:.tex=.figures)
.PHONY : figures

%.dvi : %.tex %.figures
	if ( latex $< ) ; then latex $< ; \
	else perl -e '($$rt,$$wt)=(stat"$*.tex")[8,9];utime($$rt,$$wt,"$*.dvi")' ; echo ; fi

# parse a .tex file and extract the .eps files
%.figures : %.tex
	perl -e 'print"$*.figures : ";while(<>){print "$$+ "if s/^[^%]*\{([^[\}|\{|%]*\.eps)\}/\1/}print"\n"' $< > $@

%.eps : %.agr
	gracebat -printfile $@ $<
	perl -pi -e 's/%%BoundingBox:.*/%%BoundingBox: $(GraceBoundingBox)/' $@

%.eps : %.ps
	ps2epsi $< $@

%.ps : %.agr
	gracebat -printfile $@ $<

%.ps : %.dvi
	dvips $< -o $@

# this needs some work - does pdflatex write to the .pdf file even on error?
%.pdf : %.tex %.figures
	if ( pdflatex $< ) ; then pdflatex $< ; \
	else perl -e '($$rt,$$wt)=(stat"$*.tex")[8,9];utime($$rt,$$wt,"$*.pdf")' ; echo ; fi

#%.pdf : %.dvi
#	dvipdf $< -o $@

%.tar : %.figures %.tex
	tar cf $@ $*.tex `cat $*.figures | sed 's/.*://'`

%.gz : %.tar
	gzip -f $<

%.tar.gz : %.figures %.tex
	tar cf $*.tar $*.tex `cat $*.figures | sed 's/.*://'`
	gzip -f $*.tar

clean :
	rm -f $(DVIFILES) $(PSFILES) $(PDFFILES) $(TEXSRC:%.tex=%.log) $(TEXSRC:%.tex=%.aux)
