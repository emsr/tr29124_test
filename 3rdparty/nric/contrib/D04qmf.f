c	$Log:	D04qmf.f,v $
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
c        Daubechies' Interval  QMF's with lf = 4                            c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c
	 parameter(d4h1 = 0.4829629131445341d0)
	 parameter(d4h2 = 0.8365163037378077d0)
	 parameter(d4h3 = 0.2241438680420134d0)
	 parameter(d4h4 =-0.1294095225512603d0)
c
	 parameter(d4g1 =  d4h4)
	 parameter(d4g2 = -d4h3)
	 parameter(d4g3 =  d4h2)
	 parameter(d4g4 = -d4h1)
c
c        left edge coefficients for the wavelets on the interval
c
	 parameter(d4h1l1 = 6.03332511928053d-01)
	 parameter(d4h1l2 = 6.90895531839104d-01)
	 parameter(d4h1l3 =-3.98312997698228d-01)
c
	 parameter(d4g1l1 =-7.96543516912183d-01)
	 parameter(d4g1l2 = 5.46392713959015d-01)
	 parameter(d4g1l3 =-2.58792248333818d-01)
c
	 parameter(d4h2l1 = 3.75174604524466d-02)
	 parameter(d4h2l2 = 4.57327659851769d-01)
	 parameter(d4h2l3 = 8.50088102549165d-01)
	 parameter(d4h2l4 = 2.23820356983114d-01)
	 parameter(d4h2l5 =-1.29222743354319d-01)
c
	 parameter(d4g2l1 =-1.00372245644139d-02)
	 parameter(d4g2l2 =-1.22351043116799d-01)
	 parameter(d4g2l3 =-2.27428111655837d-01)
	 parameter(d4g2l4 = 8.36602921223654d-01)
	 parameter(d4g2l5 =-4.83012921773304d-01)
c
c
c        right edge coefficients for the wavelets on the interval
c
	 parameter(d4h1r1 = 8.70508753349866d-01)
	 parameter(d4h1r2 = 4.34896997965703d-01)
	 parameter(d4h1r3 = 2.30389043796969d-01)
c
	 parameter(d4g1r1 = 2.57512919478482d-01)
	 parameter(d4g1r2 =-8.01422961990337d-01)
	 parameter(d4g1r3 = 5.39822500731772d-01)
c
	 parameter(d4h2r1 =-1.94233407427412d-01)
	 parameter(d4h2r2 = 1.90151418429955d-01)
	 parameter(d4h2r3 = 3.74955331645687d-01)
	 parameter(d4h2r4 = 7.67556669298114d-01)
	 parameter(d4h2r5 = 4.43149049637559d-01)
c
	 parameter(d4g2r1 = 3.71718966535296d-01)
	 parameter(d4g2r2 =-3.63906959570891d-01)
	 parameter(d4g2r3 =-7.17579999353722d-01)
	 parameter(d4g2r4 = 4.01069519430217d-01)
	 parameter(d4g2r5 = 2.31557595006790d-01)
c
c
c        left edge coefficients for the preconditioner
c
c
	 parameter(d4pc1l1 = 3.24894048898962d-01)
	 parameter(d4pc1l2 = 3.71580151158803d-02)
	 parameter(d4pc2l2 = 1.00144540498130d+00)
c
	 parameter(d4ipc1l1 = 3.07792649138669d+00)
	 parameter(d4ipc1l2 =-1.14204567242137d-01)
	 parameter(d4ipc2l2 = 9.98556681198888d-01)
c
c
c        right edge coefficients for the preconditioner
c
c
	 parameter(d4pc1r1 = 2.09629288435324d+00)
	 parameter(d4pc1r2 =-8.00813234246437d-01)
	 parameter(d4pc2r2 = 1.08984305289504d+00)
c
	 parameter(d4ipc1r1 = 4.77032578540915d-01)
	 parameter(d4ipc1r2 = 3.50522032550918d-01)
	 parameter(d4ipc2r2 = 9.17563310922261d-01)
c