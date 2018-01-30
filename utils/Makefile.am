
.DEFAULT_GOAL = utils_start

TARGETDIR=/usr/local/bin/


elprop:
	@echo ""

runelprop:
	@echo ""

utils_start:
	@echo ""


cleanall:
	find . -name \*~ -exec rm -f {} \;
	find . -name \*.*~ -exec rm -f {} \;

install:
	cp elprop.pl runelprop.pl lsttail.sh lstgrep.sh $(TARGETDIR)

uninstall:
	rm -f $(TARGETDIR)/elprop.pl $(TARGETDIR)/runelprop.pl $(TARGETDIR)/lsttail.sh $(TARGETDIR)/lstgrep.sh 
