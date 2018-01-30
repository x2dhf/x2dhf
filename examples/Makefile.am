
SUBDIRS = be bef bf bo c c+ c2 f2 gaf h2 h2+ li2 n2 ne o2 oed/h3 oed/kh oed/th+89 oed/z10-6 scmc fh dft xec

.PHONY: $(SUBDIRS)

all: $(SUBDIRS)

$(SUBDIRS):
	cd $@ && $(MAKE) clean
	cd $@ && $(MAKE) 

cleanlst:
	find . -name \*.lst -exec rm -f {} \;
	find .  -name \*~ -exec rm -f {} \;

install:
	@echo ""  	

uninstall:
	@echo ""  	

#
