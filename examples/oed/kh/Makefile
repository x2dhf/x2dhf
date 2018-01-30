
PROG = x2dhf
XHF = ../../../tests/xhf 

all: kh
	$(MAKE) clean	

kh: kh.data 
	test  -L $(PROG)  || ln -s ../../../src/$(PROG) .
	test -f $@.lst || $(XHF) $@  $@ | tee $@.lst

clean:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf $(PROG) 2> /dev/null

stop:
	$(XHF) stop

#

