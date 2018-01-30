
PROG = x2dhf
XHF = ../../../tests/xhf 

all: z10-6
	$(MAKE) clean	

z10-6: z10-6.data 
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@ $@ | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

