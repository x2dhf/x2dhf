
bin_PROGRAMS = users-guide 

# include the following for an export version of the package 
#users_guide_SOURCES = users-guide.pdf 

# otherwise use
users_guide_SOURCES = users-guide.tex users-guide 

#output_SOURCES = 

.DEFAULT_GOAL = $(guide)

texguide:
	@echo ""

pdfguide:
	@echo ""  

guide: 
	$(MAKE) pdf html

pdf:	users-guide.tex   
	pdflatex  users-guide.tex; pdflatex  users-guide.tex

html:	users-guide.tex 
	latex2html -split 0  users-guide.tex 2>/dev/null

cleanall:
	rm -f *.log *.toc *.pdf rm -rf users-guide 

install:
	@echo ""  	

uninstall:
	@echo ""  	
