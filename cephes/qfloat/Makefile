
CC = $(HOME)/bin/bin/gcc -std=c11 -DANSIPROT

INC_DIR = .
OBJ_DIR = obj

OBJS = \
  $(OBJ_DIR)/drand.o \
  $(OBJ_DIR)/qsindg.o \
  $(OBJ_DIR)/qexp2.o \
  $(OBJ_DIR)/qpsi.o \
  $(OBJ_DIR)/qfltbi.o \
  $(OBJ_DIR)/qsqrt.o \
  $(OBJ_DIR)/qplanck.o \
  $(OBJ_DIR)/qlog10.o \
  $(OBJ_DIR)/qatn.o \
  $(OBJ_DIR)/cmplx.o \
  $(OBJ_DIR)/qellpe.o \
  $(OBJ_DIR)/qlog1.o \
  $(OBJ_DIR)/qjyqn.o \
  $(OBJ_DIR)/log.o \
  $(OBJ_DIR)/qcgamma.o \
  $(OBJ_DIR)/qsin.o \
  $(OBJ_DIR)/exp.o \
  $(OBJ_DIR)/qexpn.o \
  $(OBJ_DIR)/qgamma.o \
  $(OBJ_DIR)/qincbi.o \
  $(OBJ_DIR)/qspenc.o \
  $(OBJ_DIR)/qstudt.o \
  $(OBJ_DIR)/qndtri.o \
  $(OBJ_DIR)/sin.o \
  $(OBJ_DIR)/qconst.o \
  $(OBJ_DIR)/qflti.o \
  $(OBJ_DIR)/qerf.o \
  $(OBJ_DIR)/qcosh.o \
  $(OBJ_DIR)/qparanoi.o \
  $(OBJ_DIR)/qrand.o \
  $(OBJ_DIR)/qfresf.o \
  $(OBJ_DIR)/qprob.o \
  $(OBJ_DIR)/qfac.o \
  $(OBJ_DIR)/qellie.o \
  $(OBJ_DIR)/qexp.o \
  $(OBJ_DIR)/qhyp.o \
  $(OBJ_DIR)/qndtr.o \
  $(OBJ_DIR)/qjn.o \
  $(OBJ_DIR)/qhy2f1.o \
  $(OBJ_DIR)/qexp10.o \
  $(OBJ_DIR)/qremain.o \
  $(OBJ_DIR)/qexpm1.o \
  $(OBJ_DIR)/qeuclid.o \
  $(OBJ_DIR)/qcerf.o \
  $(OBJ_DIR)/qmtst.o \
  $(OBJ_DIR)/qellpj.o \
  $(OBJ_DIR)/igam.o \
  $(OBJ_DIR)/qasin.o \
  $(OBJ_DIR)/qjypn.o \
  $(OBJ_DIR)/powi.o \
  $(OBJ_DIR)/qcpolylog.o \
  $(OBJ_DIR)/qcmplx.o \
  $(OBJ_DIR)/igami.o \
  $(OBJ_DIR)/sqrt.o \
  $(OBJ_DIR)/qbeta.o \
  $(OBJ_DIR)/qpolylog.o \
  $(OBJ_DIR)/incbi.o \
  $(OBJ_DIR)/qsimq.o \
  $(OBJ_DIR)/qigami.o \
  $(OBJ_DIR)/qshici.o \
  $(OBJ_DIR)/qcbrt.o \
  $(OBJ_DIR)/pow.o \
  $(OBJ_DIR)/floor.o \
  $(OBJ_DIR)/qcalc.o \
  $(OBJ_DIR)/qkn.o \
  $(OBJ_DIR)/incbet.o \
  $(OBJ_DIR)/qincb.o \
  $(OBJ_DIR)/qasinh.o \
  $(OBJ_DIR)/qkne.o \
  $(OBJ_DIR)/qyn.o \
  $(OBJ_DIR)/qccalc.o \
  $(OBJ_DIR)/qerfc.o \
  $(OBJ_DIR)/ndtri.o \
  $(OBJ_DIR)/const.o \
  $(OBJ_DIR)/qzetac.o \
  $(OBJ_DIR)/qremquo.o \
  $(OBJ_DIR)/qsinh.o \
  $(OBJ_DIR)/qin.o \
  $(OBJ_DIR)/qchyp1f1.o \
  $(OBJ_DIR)/qairy.o \
  $(OBJ_DIR)/mtherr.o \
  $(OBJ_DIR)/qacosh.o \
  $(OBJ_DIR)/qtan.o \
  $(OBJ_DIR)/qellik.o \
  $(OBJ_DIR)/qellpk.o \
  $(OBJ_DIR)/qei.o \
  $(OBJ_DIR)/polevl.o \
  $(OBJ_DIR)/qcos.o \
  $(OBJ_DIR)/qsici.o \
  $(OBJ_DIR)/qatanh.o \
  $(OBJ_DIR)/qpow.o \
  $(OBJ_DIR)/qlog.o \
  $(OBJ_DIR)/qdawsn.o \
  $(OBJ_DIR)/qigam.o \
  $(OBJ_DIR)/qk0.o \
  $(OBJ_DIR)/qine.o \
  $(OBJ_DIR)/gamma.o \
  $(OBJ_DIR)/qtanh.o \
  $(OBJ_DIR)/qfloor.o

libcephes_qfloat.so: $(OBJS)
	$(CC) -shared -o libcephes_qfloat.so $(OBJS)
	cp libcephes_qfloat.so libcephes_qfloat.dll
	cp libcephes_qfloat.so libcephes_qfloat.dll ..

# /([a-z0-9]*)\.o
# /\1.o: \1.c\n\t\$\(CC\) -c -fPIC -I\$\(INC_DIR\) -o \$\(OBJ_DIR\)/\1.o \1.c\n

$(OBJ_DIR)/drand.o: drand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/drand.o drand.c

$(OBJ_DIR)/qsindg.o: qsindg.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsindg.o qsindg.c

$(OBJ_DIR)/qexp2.o: qexp2.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qexp2.o qexp2.c

$(OBJ_DIR)/qpsi.o: qpsi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qpsi.o qpsi.c

$(OBJ_DIR)/qfltbi.o: qfltbi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qfltbi.o qfltbi.c

$(OBJ_DIR)/qsqrt.o: qsqrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsqrt.o qsqrt.c

$(OBJ_DIR)/qplanck.o: qplanck.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qplanck.o qplanck.c

$(OBJ_DIR)/qlog10.o: qlog10.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qlog10.o qlog10.c

$(OBJ_DIR)/qatn.o: qatn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qatn.o qatn.c

$(OBJ_DIR)/cmplx.o: cmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/cmplx.o cmplx.c

$(OBJ_DIR)/qellpe.o: qellpe.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qellpe.o qellpe.c

$(OBJ_DIR)/qlog1.o: qlog1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qlog1.o qlog1.c

$(OBJ_DIR)/qjyqn.o: qjyqn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qjyqn.o qjyqn.c

$(OBJ_DIR)/log.o: log.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/log.o log.c

$(OBJ_DIR)/qcgamma.o: qcgamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcgamma.o qcgamma.c

$(OBJ_DIR)/qsin.o: qsin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsin.o qsin.c

$(OBJ_DIR)/exp.o: exp.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/exp.o exp.c

$(OBJ_DIR)/qexpn.o: qexpn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qexpn.o qexpn.c

$(OBJ_DIR)/qgamma.o: qgamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qgamma.o qgamma.c

$(OBJ_DIR)/qincbi.o: qincbi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qincbi.o qincbi.c

$(OBJ_DIR)/qspenc.o: qspenc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qspenc.o qspenc.c

$(OBJ_DIR)/qstudt.o: qstudt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qstudt.o qstudt.c

$(OBJ_DIR)/qndtri.o: qndtri.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qndtri.o qndtri.c

$(OBJ_DIR)/sin.o: sin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sin.o sin.c

$(OBJ_DIR)/qconst.o: qconst.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qconst.o qconst.c

$(OBJ_DIR)/qflti.o: qflti.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qflti.o qflti.c

$(OBJ_DIR)/qerf.o: qerf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qerf.o qerf.c

$(OBJ_DIR)/qcosh.o: qcosh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcosh.o qcosh.c

$(OBJ_DIR)/qparanoi.o: qparanoi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qparanoi.o qparanoi.c

$(OBJ_DIR)/qrand.o: qrand.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qrand.o qrand.c

$(OBJ_DIR)/qfresf.o: qfresf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qfresf.o qfresf.c

$(OBJ_DIR)/qprob.o: qprob.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qprob.o qprob.c

$(OBJ_DIR)/qfac.o: qfac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qfac.o qfac.c

$(OBJ_DIR)/qellie.o: qellie.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qellie.o qellie.c

$(OBJ_DIR)/qexp.o: qexp.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qexp.o qexp.c

$(OBJ_DIR)/qhyp.o: qhyp.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qhyp.o qhyp.c

$(OBJ_DIR)/qndtr.o: qndtr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qndtr.o qndtr.c

$(OBJ_DIR)/qjn.o: qjn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qjn.o qjn.c

$(OBJ_DIR)/qhy2f1.o: qhy2f1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qhy2f1.o qhy2f1.c

$(OBJ_DIR)/qexp10.o: qexp10.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qexp10.o qexp10.c

$(OBJ_DIR)/qremain.o: qremain.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qremain.o qremain.c

$(OBJ_DIR)/qexpm1.o: qexpm1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qexpm1.o qexpm1.c

$(OBJ_DIR)/qeuclid.o: qeuclid.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qeuclid.o qeuclid.c

$(OBJ_DIR)/qcerf.o: qcerf.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcerf.o qcerf.c

$(OBJ_DIR)/qmtst.o: qmtst.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qmtst.o qmtst.c

$(OBJ_DIR)/qellpj.o: qellpj.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qellpj.o qellpj.c

$(OBJ_DIR)/igam.o: igam.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igam.o igam.c

$(OBJ_DIR)/qasin.o: qasin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qasin.o qasin.c

$(OBJ_DIR)/qjypn.o: qjypn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qjypn.o qjypn.c

$(OBJ_DIR)/powi.o: powi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/powi.o powi.c

$(OBJ_DIR)/qcpolylog.o: qcpolylog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcpolylog.o qcpolylog.c

$(OBJ_DIR)/qcmplx.o: qcmplx.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcmplx.o qcmplx.c

$(OBJ_DIR)/igami.o: igami.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/igami.o igami.c

$(OBJ_DIR)/sqrt.o: sqrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/sqrt.o sqrt.c

$(OBJ_DIR)/qbeta.o: qbeta.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qbeta.o qbeta.c

$(OBJ_DIR)/qpolylog.o: qpolylog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qpolylog.o qpolylog.c

$(OBJ_DIR)/incbi.o: incbi.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbi.o incbi.c

$(OBJ_DIR)/qsimq.o: qsimq.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsimq.o qsimq.c

$(OBJ_DIR)/qigami.o: qigami.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qigami.o qigami.c

$(OBJ_DIR)/qshici.o: qshici.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qshici.o qshici.c

$(OBJ_DIR)/qcbrt.o: qcbrt.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcbrt.o qcbrt.c

$(OBJ_DIR)/pow.o: pow.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/pow.o pow.c

$(OBJ_DIR)/floor.o: floor.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/floor.o floor.c

$(OBJ_DIR)/qcalc.o: qcalc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcalc.o qcalc.c

$(OBJ_DIR)/qkn.o: qkn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qkn.o qkn.c

$(OBJ_DIR)/incbet.o: incbet.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/incbet.o incbet.c

$(OBJ_DIR)/qincb.o: qincb.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qincb.o qincb.c

$(OBJ_DIR)/qasinh.o: qasinh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qasinh.o qasinh.c

$(OBJ_DIR)/qkne.o: qkne.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qkne.o qkne.c

$(OBJ_DIR)/qyn.o: qyn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qyn.o qyn.c

$(OBJ_DIR)/qccalc.o: qccalc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qccalc.o qccalc.c

$(OBJ_DIR)/qerfc.o: qerfc.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qerfc.o qerfc.c

$(OBJ_DIR)/ndtri.o: ndtri.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/ndtri.o ndtri.c

$(OBJ_DIR)/const.o: const.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/const.o const.c

$(OBJ_DIR)/qzetac.o: qzetac.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qzetac.o qzetac.c

$(OBJ_DIR)/qremquo.o: qremquo.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qremquo.o qremquo.c

$(OBJ_DIR)/qsinh.o: qsinh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsinh.o qsinh.c

$(OBJ_DIR)/qin.o: qin.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qin.o qin.c

$(OBJ_DIR)/qchyp1f1.o: qchyp1f1.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qchyp1f1.o qchyp1f1.c

$(OBJ_DIR)/qairy.o: qairy.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qairy.o qairy.c

$(OBJ_DIR)/mtherr.o: mtherr.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/mtherr.o mtherr.c

$(OBJ_DIR)/qacosh.o: qacosh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qacosh.o qacosh.c

$(OBJ_DIR)/qtan.o: qtan.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qtan.o qtan.c

$(OBJ_DIR)/qellik.o: qellik.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qellik.o qellik.c

$(OBJ_DIR)/qellpk.o: qellpk.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qellpk.o qellpk.c

$(OBJ_DIR)/qei.o: qei.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qei.o qei.c

$(OBJ_DIR)/polevl.o: polevl.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/polevl.o polevl.c

$(OBJ_DIR)/qcos.o: qcos.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qcos.o qcos.c

$(OBJ_DIR)/qsici.o: qsici.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qsici.o qsici.c

$(OBJ_DIR)/qatanh.o: qatanh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qatanh.o qatanh.c

$(OBJ_DIR)/qpow.o: qpow.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qpow.o qpow.c

$(OBJ_DIR)/qlog.o: qlog.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qlog.o qlog.c

$(OBJ_DIR)/qdawsn.o: qdawsn.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qdawsn.o qdawsn.c

$(OBJ_DIR)/qigam.o: qigam.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qigam.o qigam.c

$(OBJ_DIR)/qk0.o: qk0.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qk0.o qk0.c

$(OBJ_DIR)/qine.o: qine.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qine.o qine.c

$(OBJ_DIR)/gamma.o: gamma.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/gamma.o gamma.c

$(OBJ_DIR)/qtanh.o: qtanh.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qtanh.o qtanh.c

$(OBJ_DIR)/qfloor.o: qfloor.c
	$(CC) -c -fPIC -I$(INC_DIR) -o $(OBJ_DIR)/qfloor.o qfloor.c


clean:
	rm -f $(OBJS)
