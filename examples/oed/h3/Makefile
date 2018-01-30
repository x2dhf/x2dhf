
PROG = x2dhf
XHF = ../../../tests/xhf 

all: h3-1 h3-2
	$(MAKE) clean	

h3-1: h3-1.data 
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@   $@        | tee $@.lst

h3-2: h3-2.data 
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ h3-1 $@ | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

