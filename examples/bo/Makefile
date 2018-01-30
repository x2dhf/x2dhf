
PROG = x2dhf
XHF = ../../tests/xhf 

all: bo

bo: bo.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || (cp $@.pun gaussian.pun; cp $@.out gaussian.out; $(XHF) $@ $@-1 | tee $@.lst)
	$(MAKE) clean

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) gaussian.* 2> /dev/null

stop:
	$(XHF) stop

#

