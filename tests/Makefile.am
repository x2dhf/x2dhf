
bin_PROGRAMS = xhf

xhf_SOURCES = xhf


.DEFAULT_GOAL = help

XHF = ./xhf
EXAMPLES_DIR = ../examples


clean:
	$(MAKE) cleanbin cleandata cleanlst

cleanbin:
	rm -f *.dat *.orb *.coul *.exch stop_x2dhf

cleandata:
	rm -f *.data *.pun *.out

cleanlst:
	rm -rf *.lst

help:
	@echo ""
	@echo "Usage: make case"
	@echo "       to perform test calculations for a specific system found in ./examples directory."
	@echo ""
	@echo "       Supported values of case are: "
	@echo "         be bef bf_g94 bo c c1d2 c3p1 c+ c2a c2b c2c dft f2 "
	@echo "         fh-0 fh-m1 fh-m2 fh-p1 fh-p2 gaf h2 h2+ li2 "
	@echo "         n2 ne o2 oed-h3 oed-kh oed-th+89 oed-z10-6 scmc xec"
	@echo ""
	@echo "Type 'make all' to run all the cases."
	@echo ""

	bash -c " [[ -L x2dhf ]] ||  ln -s ../tests/x2dhf x2dhf "

stop:
	$(XHF) stop

all: be bef bf_g94 bo c c1d2 c3p1 c+ c2a c2b c2c dft f2 fh-0 fh-m1 fh-m2 fh-p1 fh-p2 gaf h2 h2+ li2 n2 ne o2 oed-h3 oed-kh oed-th+89 oed-z10-6 scmc

fall: be bef bf_g94 bo c c1d2 c3p1 c+ c2a c2b c2c dft f2 h2 li2 n2 ne o2 oed-h3 oed-kh oed-th+89 oed-z10-6 scmc

be: ../examples/be/be.data ../examples/be/be-1.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@   $@-1        | tee $@.lst
	$(XHF) $(EXAMPLES_DIR)/$@/$@-1 $@-1 $@-2 | tee $@-1.lst

bef: ../examples/bef/bef.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

bf_g94: ../examples/bf/bf_g94.data
	cp $(EXAMPLES_DIR)/bf/$@.pun gaussian.pun
	cp $(EXAMPLES_DIR)/bf/$@.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/bf/$@ $@-1 | tee $@.lst

bo: ../examples/bo/bo.data
	cp $(EXAMPLES_DIR)/bo/$@.pun gaussian.pun
	cp $(EXAMPLES_DIR)/bo/$@.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

c: ../examples/c/c.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

c1d2: ../examples/c/c1d2.data
	$(XHF) $(EXAMPLES_DIR)/c/$@ $@-1 | tee $@.lst

c3p1: ../examples/c/c1d2.data
	$(XHF) $(EXAMPLES_DIR)/c/$@ $@-1 | tee $@.lst

c+: ../examples/c+/c+.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

c2a: ../examples/c2/c2a.data
	cp $(EXAMPLES_DIR)/c2/c2.pun gaussian.pun
	cp $(EXAMPLES_DIR)/c2/c2.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/c2/$@ $@-1 | tee $@.lst

c2b: ../examples/c2/c2b.data
	$(XHF) $(EXAMPLES_DIR)/c2/$@ c2a-1 $@-1 | tee $@.lst

c2c: ../examples/c2/c2c.data
	$(XHF) $(EXAMPLES_DIR)/c2/$@ c2b-1 $@-1 | tee $@.lst


dft: ../examples/dft/be-1.data ../examples/dft/be-2.data
	$(XHF) $(EXAMPLES_DIR)/$@/be-1 be-b88  | tee be-1.lst
	$(XHF) $(EXAMPLES_DIR)/$@/be-2 be-b88 be-b88+lyp | tee be-2.lst

f2: ../examples/f2/f2.data
	cp $(EXAMPLES_DIR)/$@/$@.pun gaussian.pun
	cp $(EXAMPLES_DIR)/$@/$@.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

fh-0: ../examples/fh/fh-0.data
	$(XHF) $(EXAMPLES_DIR)/fh/$@ $@-1 | tee $@.lst

fh-m1: ../examples/fh/fh-m1.data
	$(XHF) $(EXAMPLES_DIR)/fh/$@ fh-0-1  $@-1 | tee $@.lst

fh-m2: ../examples/fh/fh-m2.data
	$(XHF) $(EXAMPLES_DIR)/fh/$@ fh-m1-1 $@-1 | tee $@.lst

fh-p1: ../examples/fh/fh-p1.data
	$(XHF) $(EXAMPLES_DIR)/fh/$@ fh-0-1  $@-1 | tee $@.lst

fh-p2: ../examples/fh/fh-p2.data
	$(XHF) $(EXAMPLES_DIR)/fh/$@ fh-p1-1 $@-1 | tee $@.lst

gaf: ../examples/gaf/gaf.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

h2: ../examples/h2/h2a.data ../examples/h2/h2b.data
	$(XHF) $(EXAMPLES_DIR)/h2/h2a h2a | tee h2a.lst
	$(XHF) $(EXAMPLES_DIR)/h2/h2b h2b | tee h2b.lst

h2+: ../examples/h2+/h2+.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

li2: ../examples/li2/li2.data
	cp $(EXAMPLES_DIR)/$@/$@.pun gaussian.pun
	cp $(EXAMPLES_DIR)/$@/$@.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

n2: ../examples/n2/n2.data
	$(XHF) $(EXAMPLES_DIR)/$@/$@ $@-1 | tee $@.lst

ne: ../examples/ne/ne-a.data ../examples/ne/ne-b.data ../examples/ne/ne-c.data
	$(XHF) $(EXAMPLES_DIR)/ne/ne-a ne-a | tee ne-a.lst
	$(XHF) $(EXAMPLES_DIR)/ne/ne-b ne-a ne-b | tee ne-b.lst
	$(XHF) $(EXAMPLES_DIR)/ne/ne-c ne-c | tee ne-c.lst

o2: ../examples/o2/o2-a.data ../examples/o2/o2-b.data
	cp $(EXAMPLES_DIR)/$@/$@.pun gaussian.pun
	cp $(EXAMPLES_DIR)/$@/$@.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/o2/o2-a o2-a | tee o2-a.lst
	$(XHF) $(EXAMPLES_DIR)/o2/o2-b o2-b | tee o2-b.lst

oed-h3: ../examples/oed/h3/h3-1.data  ../examples/oed/h3/h3-2.data
	$(XHF) $(EXAMPLES_DIR)/oed/h3/h3-1 h3-1      | tee h3-1.lst
	$(XHF) $(EXAMPLES_DIR)/oed/h3/h3-2 h3-1 h3-2 | tee h3-2.lst

oed-kh: ../examples/oed/kh/kh.data
	$(XHF) $(EXAMPLES_DIR)/oed/kh/kh kh-1  | tee kh.lst


oed-th+89: ../examples/oed/th+89/th+89_1s.data ../examples/oed/th+89/th+89_2s.data
	$(XHF) $(EXAMPLES_DIR)/oed/th+89/th+89_1s th+89_1s-1            | tee th+89_1s.lst
	$(XHF) $(EXAMPLES_DIR)/oed/th+89/th+89_2s th+89_1s-1 th+89_2s-1 | tee th+89_2s.lst

oed-z10-6: ../examples/oed/z10-6/z10-6.data
	$(XHF) $(EXAMPLES_DIR)/oed/z10-6/z10-6 z10-6-1 | tee z10-6.lst

scmc: ../examples/scmc/be-scmc.data
	$(XHF) $(EXAMPLES_DIR)/scmc/be-scmc be-scmc-1  | tee be-scmc.lst

xec: ../examples/xec/xe-c-1.data
	cp $(EXAMPLES_DIR)/xec/xe-c.pun gaussian.pun
	cp $(EXAMPLES_DIR)/xec/xe-c.out gaussian.out
	$(XHF) $(EXAMPLES_DIR)/xec/xe-c-1 xe-c-1 | tee xe-c-1.lst
	$(XHF) $(EXAMPLES_DIR)/xec/xe-c-2 xe-c-1 xe-c-2 | tee xe-c-2.lst
	$(XHF) $(EXAMPLES_DIR)/xec/xe-c-3 xe-c-2 xe-c-3 | tee xe-c-3.lst


install:
	@echo ""
