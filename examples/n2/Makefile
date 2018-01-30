
PROG = x2dhf
XHF = ../../tests/xhf 

all: n2
	$(MAKE) clean

n2: n2.data
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@-1 | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

