# Main program

add_executable(x2dhf blas.F90 card.f90 checkOrbSym.f90 checkOrtho.f90
checkPot.f90 checkSym.f90 commons16.f90 commons8.f90 contrib.f90
contriborbDFT.f90 contriborb.f90 coulAsymptCont.f90 coulAsympt.f90
coulMCSOR.f90 coulMom.f90 coulSOR.f90 dftauxfun.f90 dftex.f90 dgamit.f
diffmu.f90 diffnu.f90 discret.f90 dointerp.f90 dointerp_mu.f90
dointerp_nu.f90 doSCF.f90 doSOR.f90 Eab1DFT.f90 Eab1HF.f90 Eab2DFT.f90
EabDFT.f90 Eab.f90 EaDFT.f90 Ea.f90 eclyptot.f90 ecvwntot.f90
ehemsic.f90 es.f90 etotalDFT.f90 etotal.f90 etotalGauss.f90
etotalOrb.f90 evalGauss.f90 exbe88.f90 exbe88sup.f90
exchAsymptCont.f90 exchAsympt.f90 exchMCSOR.f90 exchMom.f90
exchSOR.f90 excont.f90 exint.f90 exocc.f90 expw86.f90 expw86sup.f90
expw91.f90 expw91sup.f90 extinorg.f90 exxalpha.f90 factor2.f90
factor.f90 fbe88.f90 fbe88sup.f90 fdften.f90 fdftpot.f90 fermi.f90
flp.f90 flypcs.f90 flypsupcs.f90 fockDFT.f90 fock.f90 fockform.f90
fpw86.f90 fpw86sup.f90 fvwncs.f90 fvwnsupcs.f90 getCpuTime.f90
getDateTime.f90 gmtran.f90 grid_old.f90 hydrogenOrbA.f90
hydrogenOrbB.f90 hypg1f1.f90 hypg2f1.f90 inCard.f90 inFloat.f90
inInt.f90 initAddr.f90 initArrays.f90 initCBlocks.f90 initDisk.f90
initExWeights.f90 initGauss.f90 initHF.f90 initHyd.f90 initMesh.f90
initOrbPot.f90 initPot.f90 initSuppl.f90 inputData.f90 inStr.f90
interpolq.f90 kh.f90 lagrpolq.f90 leexact1.f90 locenergy.f90
lpcoeff.f90 lpcoeffq.f90 lpderq.f90 mcsor.f90 memory.f90 mesh.f90
meshsize.f90 mpoleMom.f90 mulex.f90 multf4.f90 n2f.f90 nfng.f90
nmucalc.f90 nnucalc.f90 norm94.f90 norm.f90 nuclder1.f90 nuclder.f90
orbAsymptDet.f90 orbAsymptSet.f90 orbMCSOR.f90 orbSOR.f90 ortho.f90
params.f90 plaguer.f90 plegendg.f90 pmtx.f90 potAsympt.f90 pottf.f90
prepdiff1.f90 prepGauss.f90 prepSCF.f90 prepwexch.f90 printBanner.f90
printCase.f90 printElConf.f90 printResults.f90 propet2.f90 propet3.f90
prtdenA.f90 prtdenB.f90 prtmatcw.f90 prtmat.f90 prttoten.f90
putin2.f90 putin3.f90 putin4.f90 putin.f90 putout34.f90 putout.f90
qsortf.f radialden.f90 reada16.f90 reada8.f90 rexponents.f90
rfdexch.f90 rfunaux16.f90 rfunaux.f90 rfun.f90 rfun_int.f90
rheader.f90 rrec.f90 scf.f90 schedSOR.f90 scheduler.f90 scmc.f90
separator.f90 setCi.f90 setCiOrb.f90 setDefaults.f90 setmaxunit.f90
setOmega.f90 setOmegaOrb.f90 setOmegaPot2.f90 setOmegaPot.f90
setPrecision.f90 slaterp.f90 solver.f90 sor.f90 sort.f90 swap.f90
tail.f90 tden.f90 testn2f.f90 testnfng.f90 util.f90 vcoul.f90
vexch.f90 vlpcoeff.f90 vpoly1q.f90 vpolyq.f90 wrec.f90 wrecf.f90
writea128.f90 writea128f.f90 writea32.f90 writea32f.f90 writea64.f90
writea64f.f90 writea.f90 writeDisk4dd.f90 writeDisk.f90 wtdexch1.f90
wtdexch1f.f90 wtdexch.f90 wtdisk128.f90 wtdisk128f.f90 wtdisk32.f90
wtdisk32f.f90 wtdisk64.f90 wtdisk64f.f90 wtdisknat.f90 x2dhf.f90
zeroArray.f90 zeroArrays.f90 zgsz1.f90 zgsz1g.f90 zgsz2.f90 zgsz2g.f90
zz1.f90 zz1g.f90 zz2.f90 zz2g.f90 )

# Install
install (TARGETS x2dhf DESTINATION bin)
