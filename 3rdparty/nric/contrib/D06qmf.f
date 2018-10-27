c	$Log:	D06qmf.f,v $
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c        Written by M.E. Brewster and G. Beylkin                            c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c                                                                           c
c        This file is a part of                                             c
c        Double Precision Fast Wavelet Transform Library                    c
c        Contains proprietary information supplied by GB Consulting.        c
c        Copyright (C), 1992 GB Consulting. All rights reserved             c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c                                                                           c
c        Daubechies' Interval QMF's with lf = 6                             c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c
         parameter(d6h1 = 0.3326705529500825d0)
         parameter(d6h2 = 0.8068915093110924d0)
         parameter(d6h3 = 0.4598775021184914d0)
         parameter(d6h4 =-0.1350110200102546d0)
         parameter(d6h5 =-0.0854412738820267d0)
         parameter(d6h6 = 0.0352262918857095d0)
c
         parameter(d6g1 = d6h6)
         parameter(d6g2 =-d6h5)
         parameter(d6g3 = d6h4)
         parameter(d6g4 =-d6h3)
         parameter(d6g5 = d6h2)
         parameter(d6g6 =-d6h1)
c
c        left edge coefficients for the wavelets on the interval
c
	 parameter(d6h1l1 = 3.88899763761418d-01)
	 parameter(d6h1l2 =-8.82078282004781d-02)
	 parameter(d6h1l3 =-8.47841308471311d-01)
	 parameter(d6h1l4 = 3.49487436741469d-01)
c
	 parameter(d6g1l1 = 5.83780969578456d-01)
	 parameter(d6g1l2 = 7.93618840200530d-01)
	 parameter(d6g1l3 = 1.60955177164887d-01)
	 parameter(d6g1l4 =-5.88417112312549d-02)
c
	 parameter(d6h2l1 =-6.21148317851952d-01)
	 parameter(d6h2l2 = 5.22527393204687d-01)
	 parameter(d6h2l3 =-2.00008003001907d-01)
	 parameter(d6h2l4 = 3.37867348635904d-01)
	 parameter(d6h2l5 =-3.99770770411862d-01)
	 parameter(d6h2l6 = 1.64820129734343d-01)
c
	 parameter(d6g2l1 =-3.49340182426010d-01)
	 parameter(d6g2l2 = 2.98920557206519d-01)
	 parameter(d6g2l3 =-3.28301343932461d-01)
	 parameter(d6g2l4 =-3.32263713173123d-01)
	 parameter(d6g2l5 = 6.98249720023120d-01)
	 parameter(d6g2l6 =-2.87878999564209d-01)
c
	 parameter(d6h3l1 =-9.58786377870484d-03)
	 parameter(d6h3l2 = 3.71225516842464d-04)
	 parameter(d6h3l3 = 3.26009710110353d-01)
	 parameter(d6h3l4 = 8.01648164413686d-01)
	 parameter(d6h3l5 = 4.72055262027448d-01)
	 parameter(d6h3l6 =-1.40042080963951d-01)
	 parameter(d6h3l7 =-8.54251000994755d-02)
	 parameter(d6h3l8 = 3.52196236519723d-02)
c
	 parameter(d6g3l1 =-1.01505899532729d-03)
	 parameter(d6g3l2 = 3.93013301879430d-05)
	 parameter(d6g3l3 = 3.45143711309844d-02)
	 parameter(d6g3l4 = 8.48698103307436d-02)
	 parameter(d6g3l5 =-1.33730702398544d-01)
	 parameter(d6g3l6 =-4.60406454041029d-01)
	 parameter(d6g3l7 = 8.06893221779630d-01)
	 parameter(d6g3l8 =-3.32671258977905d-01)
c
c
c        right edge coefficients for the wavelets on the interval
c
	 parameter(d6h1r1 = 9.09684994311124d-01)
	 parameter(d6h1r2 = 3.82360655862306d-01)
	 parameter(d6h1r3 = 1.50987215328402d-01)
	 parameter(d6h1r4 = 5.89610106858016d-02)
c
	 parameter(d6g1r1 =-7.22194876326683d-02)
	 parameter(d6g1r2 = 4.26562218788225d-01)
	 parameter(d6g1r3 =-8.04233140972164d-01)
	 parameter(d6g1r4 = 4.07477697635819d-01)
c
	 parameter(d6h2r1 =-2.90407851090693d-01)
	 parameter(d6h2r2 = 4.18999228996535d-01)
	 parameter(d6h2r3 = 4.96964372120369d-01)
	 parameter(d6h2r4 = 4.90757830674591d-01)
	 parameter(d6h2r5 = 4.64362767365511d-01)
	 parameter(d6h2r6 = 1.91450544225960d-01)
c
	 parameter(d6g2r1 =-1.53505230705603d-01)
	 parameter(d6g2r2 = 5.22394237872772d-01)
	 parameter(d6g2r3 =-9.81980008912701d-02)
	 parameter(d6g2r4 =-7.67879574290977d-01)
	 parameter(d6g2r5 = 2.98515239695673d-01)
	 parameter(d6g2r6 = 1.23073831745202d-01)
c
	 parameter(d6h3r1 = 8.18354184018805d-02)
	 parameter(d6h3r2 =-1.58758215582660d-01)
	 parameter(d6h3r3 =-9.12473562312011d-02)
	 parameter(d6h3r4 = 6.04255818640831d-04)
	 parameter(d6h3r5 = 7.70293360943509d-02)
	 parameter(d6h3r6 = 5.20060177812491d-01)
	 parameter(d6h3r7 = 7.64259199273409d-01)
	 parameter(d6h3r8 = 3.15093823005453d-01)
c
	 parameter(d6g3r1 = 2.29477548350898d-01)
	 parameter(d6g3r2 =-4.45179444352116d-01)
	 parameter(d6g3r3 =-2.55869891183398d-01)
	 parameter(d6g3r4 = 1.69441479675110d-03)
	 parameter(d6g3r5 = 7.59876151023844d-01)
	 parameter(d6g3r6 = 1.39150326725670d-01)
	 parameter(d6g3r7 =-2.72547235184810d-01)
	 parameter(d6g3r8 =-1.12367571585129d-01)
c
c
c        left edge coefficients for the preconditioner
c
c
	 parameter(d6pc1l1 = 1.00794157907793d-01)
	 parameter(d6pc1l2 =-5.92931024415852d-01)
	 parameter(d6pc2l2 = 2.13725665320657d-01)
	 parameter(d6pc1l3 =-1.50964976106448d-02)
	 parameter(d6pc2l3 = 3.06854237778467d-02)
	 parameter(d6pc3l3 = 1.00018933290722d+00)
c
	 parameter(d6ipc1l1 = 9.92120992681742d+00)
	 parameter(d6ipc1l2 = 2.75240372115665d+01)
	 parameter(d6ipc2l2 = 4.67889524872776d+00)
	 parameter(d6ipc1l3 =-6.94679698232383d-01)
	 parameter(d6ipc2l3 =-1.43546705404309d-01)
	 parameter(d6ipc3l3 = 9.99810702932945d-01)
c
c
c        right edge coefficients for the preconditioner
c
c
	 parameter(d6pc1r1 = 5.64221252452501d+00)
	 parameter(d6pc1r2 =-6.66336794227380d+00)
	 parameter(d6pc2r2 = 1.73763179569470d+00)
	 parameter(d6pc1r3 = 2.41778315064289d+00)
	 parameter(d6pc2r3 =-4.65879830148413d-01)
	 parameter(d6pc3r3 = 1.05578252781022d+00)
c
	 parameter(d6ipc1r1 = 1.77235436569129d-01)
	 parameter(d6ipc1r2 = 6.79652000611259d-01)
	 parameter(d6ipc2r2 = 5.75495914886965d-01)
	 parameter(d6ipc1r3 =-1.05969449845818d-01)
	 parameter(d6ipc2r3 = 2.53946179271153d-01)
	 parameter(d6ipc3r3 = 9.47164755675664d-01)
c
