
PROG = x2dhf
XHF = ../../tests/xhf 

all: be be-1
	$(MAKE) clean	

be: be.data 
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@   $@-1        | tee $@.lst

be-1: be-1.data 
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@ $@-2 | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

