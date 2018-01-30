
PROG = x2dhf
XHF = ../../tests/xhf 

all: be-scmc
	$(MAKE) clean	

be-scmc: be-scmc.data 
	test  -L $(PROG)  || ln -s ../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@   $@-1        | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

