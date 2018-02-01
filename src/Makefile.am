bin_PROGRAMS = x2dhf

x2dhf_SOURCES = x2dhf.f90
.DEFAULT_GOAL = prog

TARGETDIR=/usr/local/bin/

# Modify the following three lines, if necessary
MOD_PATH =
##FFLAGS   = -O -g -fcheck=all $(MOD_PATH)
##FFLAGS   = -O -g  -fbounds-check  $(MOD_PATH)
##FFLAGS   = -O2  -fforce-addr -fstrength-reduce -funroll-loops
##FFLAGS   = -O2   -ffpe-trap=underflow
FC       = gfortran
FFLAGS   = -O2 -Wall

AM_FLAGS   =

#LIBS =-static-libgfortran
LIBS =

MODS = params.o discret.o card.o commons8.o commons16.o kh.o memory.o scf.o scheduler.o solver.o

OBJS = x2dhf.o \
     blas.o \
     checkOrbSym.o \
     checkOrtho.o \
     checkPot.o \
     checkSym.o \
     contrib.o \
     contriborb.o \
     contriborbDFT.o \
     coulAsympt.o \
     coulAsymptCont.o \
     coulMCSOR.o \
     coulMom.o \
     coulSOR.o \
     dftex.o \
     dftauxfun.o \
     dgamit.o \
     diffmu.o \
     diffnu.o \
     dointerp.o \
     dointerp_mu.o \
     dointerp_nu.o \
     doSCF.o \
     doSOR.o \
     Ea.o \
     EaDFT.o \
     Eab.o \
     Eab1HF.o \
     Eab1DFT.o \
     Eab2DFT.o \
     EabDFT.o \
     eclyptot.o \
     ecvwntot.o \
     ehemsic.o \
     es.o \
     etotal.o \
     etotalGauss.o \
     etotalOrb.o \
     etotalDFT.o \
     evalGauss.o \
     exbe88.o \
     exbe88sup.o \
     exchAsympt.o \
     exchAsymptCont.o \
     exchMCSOR.o \
     exchMom.o \
     exchSOR.o \
     excont.o \
     exint.o \
     exocc.o \
     expw86.o \
     expw86sup.o \
     expw91.o \
     expw91sup.o \
     extinorg.o \
     exxalpha.o \
     factor.o \
     factor2.o \
     fbe88.o \
     fbe88sup.o \
     fdften.o \
     fdftpot.o \
     fermi.o \
     flp.o \
     flypcs.o \
     flypsupcs.o \
     fock.o \
     fockDFT.o \
     fockform.o \
     fpw86.o \
     fpw86sup.o \
     fvwncs.o \
     fvwnsupcs.o \
     getCpuTime.o \
     getDateTime.o \
     grid_old.o \
     gmtran.o \
     hydrogenOrbA.o \
     hydrogenOrbB.o \
     hypg1f1.o \
     hypg2f1.o \
     inCard.o \
     inFloat.o \
     inInt.o \
     inStr.o \
     initAddr.o \
     initArrays.o \
     initCBlocks.o \
     initDisk.o \
     initExWeights.o \
     initGauss.o \
     initHF.o \
     initHyd.o \
     initMesh.o \
     initOrbPot.o \
     initPot.o \
     initSuppl.o \
     inputData.o \
     interpolq.o \
     lagrpolq.o \
     leexact1.o \
     locenergy.o \
     lpcoeff.o \
     lpcoeffq.o \
     lpderq.o \
     mcsor.o \
     mesh.o \
     meshsize.o \
     mpoleMom.o \
     mulex.o \
     multf4.o \
     n2f.o \
     nfng.o \
     nmucalc.o \
     nnucalc.o \
     norm.o \
     norm94.o \
     nuclder.o \
     nuclder1.o \
     orbAsymptDet.o \
     orbAsymptSet.o \
     orbMCSOR.o \
     orbSOR.o \
     ortho.o \
     plaguer.o \
     plegendg.o \
     pmtx.o \
     pottf.o \
     potAsympt.o \
     prepdiff1.o \
     prepGauss.o \
     prepGaussCust.o \
     prepSCF.o \
     prepwexch.o \
     printBanner.o \
     printCase.o \
     printElConf.o \
     printResults.o \
     propet2.o \
     propet3.o \
     prtdenA.o \
     prtdenB.o \
     prtmat.o \
     prtmatcw.o \
     prttoten.o \
     putin.o \
     putin2.o \
     putin3.o \
     putin4.o \
     putout.o \
     putout34.o \
     qsortf.o \
     radialden.o \
     reada8.o \
     reada16.o \
     rexponents.o \
     rfdexch.o \
     rfun.o \
     rfunaux.o \
     rfunaux16.o \
     rfun_int.o \
     rheader.o \
     rrec.o \
     schedSOR.o \
     scmc.o \
     separator.o \
     setCi.o \
     setCiOrb.o \
     setDefaults.o \
     setOmega.o \
     setOmegaOrb.o \
     setOmegaPot.o \
     setOmegaPot2.o \
     setPrecision.o \
     setmaxunit.o \
     slaterp.o \
     sor.o \
     tail.o \
     tden.o \
     testn2f.o \
     testnfng.o \
     util.o \
     vcoul.o \
     vexch.o \
     vlpcoeff.o \
     vpolyq.o \
     vpoly1q.o \
     wrec.o \
     wrecf.o \
     writea.o \
     writea32.o \
     writea32f.o \
     writea64.o \
     writea64f.o \
     writea128.o \
     writea128f.o \
     writeDisk.o \
     writeDisk4dd.o \
     wtdexch.o \
     wtdexch1.o \
     wtdexch1f.o \
     wtdisk32.o \
     wtdisk32f.o \
     wtdisk64.o \
     wtdisk64f.o \
     wtdisk128.o \
     wtdisk128f.o \
     wtdisknat.o \
     zeroArray.o \
     zeroArrays.o \
     zz1.o \
     zz2.o \
     zz1g.o \
     zz2g.o \
     zgsz1.o\
     zgsz2.o \
     zgsz1g.o \
     zgsz2g.o

prog : $(MODS) $(OBJS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS) $(LIBS) -o x2dhf



help:
	@echo ""
	@echo "Usage: make [help]"


	@echo "make"
	@echo "     compiles the program using the default compiler comp, i.e. the one set by running"
	@echo "     ./configure --enable-comp=yes"
	@echo ""

params.f90.o : params.f90
discret.f90.o : discret.f90   params.mod
card.f90.o : card.f90   params.mod
commons8.f90.o : commons8.f90   params.mod
commons16.f90.o : commons16.f90   params.mod
memory.f90.o : memory.f90
scf.f90.o : scf.f90   params.mod
scheduler.f90.o : scheduler.f90   params.mod
solver.f90.o : solver.f90   params.mod
blas.f90.o : blas.f90   params.mod
checkOrbSym.f90.o : checkOrbSym.f90   params.mod  commons8.mod
checkOrtho.f90.o : checkOrtho.f90   params.mod  commons8.mod
checkOrb.f90.o : checkOrb.f90   params.mod  discret.mod commons8.mod
checkPot.f90.o : checkPot.f90   params.mod  commons8.mod
checkSym.f90.o : checkSym.f90   params.mod  commons8.mod
coulAsympt.f90.o : coulAsympt.f90   params.mod  scf.mod  commons8.mod
coulAsymptCont.f90.o : coulAsymptCont.f90   params.mod  scf.mod  commons8.mod
coulMCSOR.f90.o : coulMCSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
coulMom.f90.o : coulMom.f90   params.mod  scf.mod  commons8.mod
contrib.f90.o : contrib.f90   params.mod  commons8.mod
contriborb.f90.o : contriborb.f90   params.mod  scf.mod  commons8.mod
contriborbDFT.f90.o : contriborbDFT.f90   params.mod  commons8.mod
coulSOR.f90.o : coulSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
dftauxfun.f90.o : dftauxfun.f90   params.mod  commons8.mod
diffmu.f90.o : diffmu.f90   params.mod  discret.mod  commons8.mod
diffnu.f90.o : diffnu.f90   params.mod  discret.mod  commons8.mod
dointerp.f90.o : dointerp.f90   params.mod  commons8.mod  commons16.mod
dointerp_mu.f90.o : dointerp_mu.f90   params.mod  commons8.mod  commons16.mod
dointerp_nu.f90.o : dointerp_nu.f90   params.mod  commons8.mod  commons16.mod
doSCF.f90.o : doSCF.f90   params.mod  commons8.mod
doSOR.f90.o : doSOR.f90   params.mod  commons8.mod
Ea.f90.o : Ea.f90   params.mod  scf.mod  commons8.mod
EaDFT.f90.o : EaDFT.f90   params.mod  memory.mod  scf.mod  commons8.mod
Eab.f90.o : Eab.f90   params.mod  scf.mod  commons8.mod
Eab1HF.f90.o : Eab1HF.f90   params.mod  scf.mod  commons8.mod
EabDFT.f90.o : EabDFT.f90   params.mod  scf.mod  commons8.mod
Eab1DFT.f90.o : Eab1DFT.f90   params.mod  memory.mod  scf.mod  commons8.mod
Eab2DFT.f90.o : Eab2DFT.f90   params.mod  scf.mod  commons8.mod
eclyptot.f90.o : eclyptot.f90   params.mod  commons8.mod
ecvwntot.f90.o : ecvwntot.f90   params.mod  commons8.mod
ehemsic.f90.o : ehemsic.f90   params.mod  commons8.mod
es.f90.o : es.f90   params.mod
etotal.f90.o : etotal.f90   params.mod  commons8.mod
etotalDFT.f90.o : etotalDFT.f90   params.mod  commons8.mod  params.mod  commons8.mod
etotalGauss.f90.o : etotalGauss.f90   params.mod  commons8.mod
evalGauss.f90.o : initGauss.f90   params.mod  commons8.mod
exbe88.f90.o : exbe88.f90   params.mod  commons8.mod
exbe88sup.f90.o : exbe88sup.f90   params.mod  commons8.mod
exchAsympt.f90.o : exchAsympt.f90   params.mod  commons8.mod
exchAsymptCont.f90.o : exchAsymptCont.f90   params.mod  commons8.mod
exchMCSOR.f90.o : exchMCSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
exchMom.f90.o : exchMom.f90   params.mod  scf.mod  commons8.mod
exchSOR.f90.o : exchSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
excont.f90.o : excont.f90   params.mod  commons8.mod
exocc.f90.o : exocc.f90   params.mod  commons8.mod
expw86.f90.o : expw86.f90   params.mod  commons8.mod
expw86sup.f90.o : expw86sup.f90   params.mod  commons8.mod
expw91sup.f90.o : expw91sup.f90   params.mod  commons8.mod
exint.f90.o : exint.f90   params.mod  commons8.mod
extinorg.o : extinorg.f90 params.mod  commons8.mod
exxalpha.f90.o : exxalpha.f90   params.mod  commons8.mod
factor.f90.o : factor.f90   params.mod
factor2.f90.o : factor2.f90   params.mod
fbe88.f90.o : fbe88.f90   params.mod  commons8.mod
fbe88sup.f90.o : fbe88sup.f90   params.mod  commons8.mod
fdften.f90.o : fdften.f90   params.mod  commons8.mod
fdftpot.f90.o : fdftpot.f90   params.mod  commons8.mod
fermi.f90.o : fermi.f90   params.mod  commons8.mod
flp.f90.o : flp.f90   params.mod
flypcs.f90.o : flypcs.f90   params.mod  commons8.mod
flypsupcs.f90.o : flypsupcs.f90   params.mod  commons8.mod
fock.f90.o : fock.f90   params.mod  scf.mod  commons8.mod
fockDFT.f90.o : fockDFT.f90   params.mod  memory.mod  scf.mod  commons8.mod
fockform.f90.o : fockform.f90   params.mod  scf.mod  commons8.mod
fpw86.f90.o : fpw86.f90   params.mod  commons8.mod
fvwncs.f90.o : fvwncs.f90   params.mod  commons8.mod
fvwnsupcs.f90.o : fvwnsupcs.f90   params.mod  commons8.mod
getCpuTime.f90.o : getCpuTime.f90   params.mod
gmtran.f90.o : gmtran.f90   params.mod
grid_old.f90.o : grid_old.f90   params.mod  discret.mod  scf.mod  commons8.mod
hydrogenOrbA.f90.o : hydrogenOrbA.f90   params.mod  commons8.mod
hydrogenOrbB.f90.o : hydrogenOrbB.f90   params.mod  commons8.mod
inCard.f90.o : inCard.f90   params.mod  card.mod  params.mod
inFloat.f90.o : inFloat.f90   params.mod  card.mod
inInt.f90.o : inInt.f90   params.mod  card.mod
initAddr.f90.o : initAddr.f90   params.mod  memory.mod  commons8.mod
initArrays.f90.o : initArrays.f90   params.mod  commons8.mod
initCBlocks.f90.o : initCBlocks.f90   params.mod  discret.mod  scf.mod  commons8.mod
initDisk.f90.o : initDisk.f90   params.mod  commons8.mod
initExWeights.f90.o : initExWeights.f90   params.mod  scf.mod  commons8.mod
initGauss.f90.o : initGauss.f90   params.mod  commons8.mod
initHF.f90.o : initHF.f90   params.mod  commons8.mod
initHyd.f90.o : initHyd.f90   params.mod  commons8.mod
initMesh.f90.o : initMesh.f90   params.mod  memory.mod  commons8.mod
initPot.f90.o : initPot.f90   params.mod  memory.mod  commons8.mod
initOrbPot.f90.o : initOrbPot.f90   params.mod  commons8.mod
initSuppl.f90.o : initSuppl.f90   params.mod  discret.mod  kh.mod  scf.mod  commons8.mod
inputData.f90.o : inputData.f90   params.mod  scf.mod  commons8.mod
inStr.f90.o : inStr.f90   params.mod  card.mod
interpolq.f90.o : interpolq.f90   params.mod  card.mod
lagrpolq.o : lagrpolq.f90  params.mod
leexact1.o : leexact1.f90 params.mod discret.mod commons8.mod
lpcoeff.f90.o : lpcoeff.f90   params.mod
lpcoeffq.f90.o : lpcoeffq.f90   params.mod
lpderq.f90.o : lpderq.f90   params.mod
locenergy.f90.o : locenergy.f90   params.mod  scf.mod  commons8.mod  params.mod  commons8.mod
mcsor.f90.o : mcsor.f90   params.mod  commons8.mod
mesh.f90.o : mesh.f90   params.mod  commons8.mod
meshsize.f90.o : meshsize.f90   params.mod  commons8.mod
mpoleMom.o : mpoleMom.f90 params.mod commons8.mod
mulex.f90.o : mulex.f90   params.mod  commons8.mod
multf4.f90.o : multf4.f90   params.mod  commons8.mod
n2f.f90.o : n2f.f90   params.mod  commons8.mod  params.mod  commons8.mod
nfng.f90.o : nfng.f90   params.mod  commons8.mod  params.mod  commons8.mod
nmucalc.f90.o : nmucalc.f90   params.mod  discret.mod  scf.mod  commons8.mod
nnucalc.f90.o : nnucalc.f90   params.mod  commons8.mod
norm.f90.o : norm.f90   params.mod  commons8.mod
norm94.f90.o : norm94.f90     params.mod  commons8.mod
nuclder.f90.o : nuclder.f90   params.mod  discret.mod  commons8.mod
nuclder1.f90.o : nuclder1.f90   params.mod  discret.mod  commons8.mod
ortho.f90.o : ortho.f90   params.mod  commons8.mod
orbAsymptDet.f90.o : orbAsymptDet.f90   params.mod  commons8.mod
orbAsymptSet.f90.o : orbAsymptSet.f90    params.mod  commons8.mod
orbMCSOR.f90.o : orbMCSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
orbSOR.f90.o : orbSOR.f90   params.mod  discret.mod  scf.mod  commons8.mod
plaguer.f90.o : plaguer.f90   params.mod
plegendg.f90.o : plegendg.f90   params.mod  commons8.mod
pmtx.f90.o : pmtx.f90   params.mod  params.mod
pottf.f90.o : pottf.f90   params.mod
potAsympt.f90.o : potAsympt.f90   params.mod  commons8.mod
prepdiff1.f90.o : prepdiff1.f90   params.mod  commons8.mod
prepGauss.f90.o : prepGauss.f90   params.mod  commons8.mod
prepGaussCust.f90.o : prepGaussCust.f90   params.mod  commons8.mod
prepSCF.f90.o : prepSCF.f90   params.mod  commons8.mod
prepwexch.f90.o : prepwexch.f90   params.mod  commons8.mod
prtdenA.f90.o : prtdenA.f90   params.mod commons8.mod
prtdenB.f90.o : prtdenB.f90   params.mod commons8.mod
prtmat.f90.o : prtmat.f90   params.mod
prtmatcw.f90.o : prtmatcw.f90   params.mod
printCase.f90.o : printCase.f90   params.mod  discret.mod  memory.mod  scf.mod  commons8.mod
printElConf.f90.o : printElConf.f90   params.mod commons8.mod
printResults.f90.o : printResults.f90   params.mod  commons8.mod
propet2.f90.o : propet2.f90   params.mod  scf.mod  commons8.mod
propet3.f90.o : propet3.f90   params.mod  scf.mod  commons8.mod
prttoten.f90.o : prttoten.f90   params.mod  commons8.mod
putin.f90.o : putin.f90   params.mod
putin2.f90.o : putin2.f90   params.mod  solver.mod
putin3.f90.o : putin3.f90   params.mod  solver.mod
putin4.f90.o : putin4.f90   params.mod  solver.mod
putout.f90.o : putout.f90   params.mod
putout34.f90.o : putout34.f90   params.mod  commons8.mod
radialden.f90.o : radialden.f90   params.mod  commons8.mod
reada8.f90.o : reada8.f90   params.mod
rexponents.f90.o : rexponents.f90   params.mod  commons8.mod
reada16.f90.o : reada16.f90   params.mod
rfdexch.f90.o : rfdexch.f90   params.mod  commons8.mod
rfun.f90.o : rfun.f90   params.mod  memory.mod  scf.mod  commons8.mod
rfunaux.f90.o : rfunaux.f90   params.mod commons8.mod
rfunaux16.f90.o : rfunaux16.f90 params.mod commons8.mod
rfun_int.f90.o : rfun_int.f90   params.mod  scf.mod  commons8.mod
rheader.f90.o : rheader.f90   params.mod  scf.mod  commons8.mod
rrec.f90.o : rrec.f90   params.mod  commons8.mod
schedSOR.f90.o : schedSOR.f90   params.mod  commons8.mod  scheduler.mod
scmc.f90.o : scmc.f90   params.mod  commons8.mod
setCi.f90.o : setCi.f90   params.mod  commons8.mod
setCiOrb.f90.o : setCiOrb.f90   params.mod  commons8.mod
setDefaults.f90.o : setDefaults.f90   params.mod  scf.mod  commons8.mod
setmaxunit.f90.o : setmaxunit.f90   commons8.mod
setOmega.f90.o : setOmega.f90   params.mod  commons8.mod
setOmegaOrb.f90.o : setOmegaOrb.f90   params.mod  commons8.mod
setOmegaPot.f90.o : setOmegaPot.f90   params.mod  commons8.mod
setOmegaPot2.f90.o : setOmegaPot2.f90   params.mod  commons8.mod
setPrecision.f90.o : setPrecision.f90   commons8.mod
slaterp.f90.o : slaterp.f90   params.mod  commons8.mod
sor.f90.o : sor.f90   params.mod  commons8.mod
tail.f90.o : tail.f90   params.mod  commons8.mod
tden.f90.o : tden.f90   params.mod  commons8.mod
testn2f.o : testn2f.f90 params.mod  commons8.mod
testnfng.o : testnfng.f90 params.mod  commons8.mod
util.f90.o : util.f90   params.mod  params.mod  params.mod  params.mod  params.mod
vcoul.f90.o : vcoul.f90   params.mod  scf.mod  commons8.mod
vexch.f90.o : vexch.f90   params.mod  scf.mod  commons8.mod
vlpcoeff.f90.o : vlpcoeff.f90   params.mod
vpolyq.f90.o : vpolyq.f90   params.mod
vpoly1q.f90.o : vpoly1q.f90   params.mod
wrec.f90.o : wrec.f90   params.mod  commons8.mod
wrecf.f90.o : wrecf.f90   params.mod  commons8.mod
writea.f90.o : writea.f90   params.mod
writea32.f90.o : writea32.f90  params.mod
writea32f.f90.o : writea32f.f90  params.mod commons8.mod
writea64.f90.o : writea64.f90  params.mod
writea64f.f90.o : writea64f.f90  params.mod commons8.mod
writea128.f90.o : writea128.f90   params.mod  commons8.mod
writea128f.f90.o : writea128f.f90   params.mod  commons8.mod
writeDisk.f90.o : writeDisk.f90   params.mod  scf.mod  commons8.mod
writeDisk4dd.f90.o : writeDisk4dd.f90   params.mod  scf.mod  commons8.mod
wtdexch.f90.o : wtdexch.f90   params.mod  commons8.mod
wtdexch1.f90.o : wtdexch1.f90   params.mod  commons8.mod
wtdexch1f.f90.o : wtdexch1f.f90   params.mod  commons8.mod
wtdisk32.f90.o : wtdisk32.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisk32f.f90.o : wtdisk32f.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisk64.f90.o : wtdisk64.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisk64f.f90.o : wtdisk64f.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisk128.f90.o : wtdisk128.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisk128f.f90.o : wtdisk128f.f90   params.mod  memory.mod  scf.mod  commons8.mod
wtdisknat.f90.o : wtdisknat.f90   params.mod  memory.mod  scf.mod  commons8.mod
zeroArray.f90.o : zeroArray.f90   params.mod
zeroArrays.f90.o : zeroArrays.f90   params.mod  memory.mod  scf.mod
zz1.f90.o : zz1.f90   params.mod  commons8.mod
zz2.f90.o : zz2.f90   params.mod  commons8.mod
zz1g.f90.o : zz1g.f90   params.mod  commons8.mod
zz2g.f90.o : zz2g.f90   params.mod  commons8.mod
zgsz1.f90.o : zgsz1.f90   params.mod  commons8.mod
zgsz2.f90.o : zgsz2.f90   params.mod  commons8.mod
zgsz1g.f90.o : zgsz1g.f90   params.mod  commons8.mod
zgsz2g.f90.o : zgsz2g.f90   params.mod  commons8.mod
mpoleMom.f90.o : mpoleMom.f90   params.mod  commons8.mod
x2dhf.f90.o : x4dhf.f90   params.mod  memory.mod  scf.mod  commons8.mod


dgamit.o : dgamit.f
	$(FC) -c $(FFLAGS) $<
qsortf.o : qsortf.f
	$(FC) -c $(FFLAGS) $<

setPrecision.o : setPrecision.f90 params.mod commons8.mod
	$(FC) -c -g -O1 $<

# rfdexch ea eab eadft eabdft etotal etotaldft prttoten



clean :
	\rm -f $(OBJS) $(MODS) *.mod *.gcov *.gcda *.gcno *.out

SUFFIXES = .f90
.f90.o:
	$(FC) -c $(FFLAGS) $<


#################################################################################

.PHONY: cleanall cleanobj cleanprog

cleanall : cleanobj cleanprog
	$(MAKE) cleanobj cleanprog cleanraw

cleanobj :
	rm -f *.o *.mod

cleanprog :
	rm -f x2dhf


install:
	cp x2dhf $(TARGETDIR)

uninstall:
	rm -f $(TARGETDIR)/x2dhf

####################################################################################################
