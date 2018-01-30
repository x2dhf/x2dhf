# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign 1.4

SUBDIRS = docs examples src tests utils 

.PHONY: $(SUBDIRS)

cleanall: $(SUBDIRS)

$(SUBDIRS):
	cd $@ && $(MAKE) cleanall
#


VERSION = f90
NAME = ./x2dhf-$(VERSION)
HOST=ameryk.fizyka.umk.pl

tgz:
#	$(MAKE) clean
	cd docs; make users-guide; cd -
	cd ..; tar  --exclude='.deps' --exclude='.svn' --exclude='src/*.o' --exclude='src/x2dhf' -czvf $(NAME)/x2dhf-$(VERSION).tgz \
	   $(NAME)/docs/Makefile.am $(NAME)/docs/users-guide.pdf $(NAME)/examples \
	   $(NAME)/src $(NAME)/tests/Makefile.am $(NAME)/tests/xhf $(NAME)/utils \
	   $(NAME)/Makefile.am $(NAME)/configure $(NAME)/install-sh $(NAME)/configure.ac $(NAME)/autoconf.sh \
	   $(NAME)/AUTHORS $(NAME)/ChangeLog $(NAME)/COPYING $(NAME)/HINTS $(NAME)/INSTALL \
	   $(NAME)/NEWS $(NAME)/README $(NAME)/TODO 


clean:
	find .  \( -name *~ -or -name \#*\# -or -name \*~ \)  -exec rm -f {} \;

upload:
	scp $(NAME).tgz  $(HOST):~/public_html/software/2dhf/
#
#
