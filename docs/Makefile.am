docfiles = users-guide.tex

if HAVE_PDFLATEX
MANNAME = users-guide
MANTEX = $(MANNAME).tex
MANAUX = $(MANNAME).aux
MANPDF = $(MANNAME).pdf
docfiles += $(MANPDF)

CLEANFILES = $(MANPDF) $(MANNAME).log $(MANNAME).idx $(MANNAME).out $(MANNAME).toc $(MANAUX)

$(MANPDF): $(srcdir)/$(MANTEX)
	$(PDFLATEX) $<
	$(PDFLATEX) $<
endif

doc_DATA = $(docfiles)
