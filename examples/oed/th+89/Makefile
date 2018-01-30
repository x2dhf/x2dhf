
PROG = x2dhf
XHF = ../../../tests/xhf 

all: th+89_1s th+89_2s
	$(MAKE) clean	

th+89_1s: th+89_1s.data
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@ | tee $@.lst

th+89_2s: th+89_2s.data
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ th+89_1s $@ | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

