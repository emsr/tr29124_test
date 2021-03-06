c	$Log:	D12qmf.f,v $
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
c        Daubechies' Interval QMF's with lf = 12                            c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c
         parameter(d12h1 =-7.8007083250323804142d-03)
         parameter(d12h2 = 1.7677118642540077410d-03)
         parameter(d12h3 = 4.4724901770781384662d-02)
         parameter(d12h4 =-2.1060292512370847991d-02)
         parameter(d12h5 =-7.2637522786376583464d-02)
         parameter(d12h6 = 3.3792942172816583271d-01)
         parameter(d12h7 = 7.8764114102865099607d-01)
         parameter(d12h8 = 4.9105594192797373304d-01)
         parameter(d12h9 =-4.8311742585698054971d-02)
         parameter(d12h10=-1.1799011114852002540d-01)
         parameter(d12h11= 3.4907120842221625153d-03)
         parameter(d12h12= 1.5404109327044824299d-02)
c
         parameter(d12g1 = d12h12)
         parameter(d12g2 =-d12h11)
         parameter(d12g3 = d12h10)
         parameter(d12g4 =-d12h9 )
         parameter(d12g5 = d12h8 )
         parameter(d12g6 =-d12h7 )
         parameter(d12g7 = d12h6 )
         parameter(d12g8 =-d12h5 )
         parameter(d12g9 = d12h4 )
         parameter(d12g10=-d12h3 )
         parameter(d12g11= d12h2 )
         parameter(d12g12=-d12h1 )
c
         parameter(d12r1 = d12g12)
         parameter(d12r2 = d12g11)
         parameter(d12r3 = d12g10)
         parameter(d12r4 = d12g9 )
         parameter(d12r5 = d12g8 )
         parameter(d12r6 = d12g7 )
         parameter(d12r7 = d12g6 )
         parameter(d12r8 = d12g5 )
         parameter(d12r9 = d12g4 )
         parameter(d12r10= d12g3 )
         parameter(d12r11= d12g2 )
         parameter(d12r12= d12g1 )
c
         parameter(d12s1 = d12h12)
         parameter(d12s2 = d12h11)
         parameter(d12s3 = d12h10)
         parameter(d12s4 = d12h9 )
         parameter(d12s5 = d12h8 )
         parameter(d12s6 = d12h7 )
         parameter(d12s7 = d12h6 )
         parameter(d12s8 = d12h5 )
         parameter(d12s9 = d12h4 )
         parameter(d12s10= d12h3 )
         parameter(d12s11= d12h2 )
         parameter(d12s12= d12h1 )
c
c        left edge coefficients for the wavelets on the interval
c
c
c
        parameter(d12h1l1 =  9.23118445982714d-01)
        parameter(d12h1l2 =  3.78106450925123d-01)
        parameter(d12h1l3 =  6.81597239670314d-02)
        parameter(d12h1l4 = -1.06339293459044d-02)
        parameter(d12h1l5 = -9.72481979806043d-03)
        parameter(d12h1l6 =  1.29788731810786d-03)
        parameter(d12h1l7 =  5.72375742644921d-03)
c
        parameter(d12h2l1 = -2.86892412498352d-01)
        parameter(d12h2l2 =  5.94046938655877d-01)
        parameter(d12h2l3 =  6.60960433566934d-01)
        parameter(d12h2l4 =  3.41298800202474d-01)
        parameter(d12h2l5 =  9.71687707534228d-02)
        parameter(d12h2l6 = -1.86085557714953d-02)
        parameter(d12h2l7 = -4.01370000538602d-02)
        parameter(d12h2l8 =  1.53730748132432d-03)
        parameter(d12h2l9 =  6.78396038981273d-03)
c
        parameter(d12h3l1 =  1.47287093229567d-01)
        parameter(d12h3l2 = -3.46894512327318d-01)
        parameter(d12h3l3 =  6.22104377284324d-02)
        parameter(d12h3l4 =  4.93802210833975d-01)
        parameter(d12h3l5 =  5.84910589484882d-01)
        parameter(d12h3l6 =  4.59878883289141d-01)
        parameter(d12h3l7 =  2.27422452903542d-01)
        parameter(d12h3l8 = -2.78895367980550d-02)
        parameter(d12h3l9 = -6.36794302334895d-02)
        parameter(d12h3l10 =  2.17771156771597d-03)
        parameter(d12h3l11 =  9.60998967044342d-03)
c
        parameter(d12h4l1 = -8.32413802141548d-02)
        parameter(d12h4l2 =  2.03563092109572d-01)
        parameter(d12h4l3 = -8.28913195517805d-02)
        parameter(d12h4l4 = -1.93990410468477d-01)
        parameter(d12h4l5 = -1.90845796716022d-02)
        parameter(d12h4l6 =  2.63640354008395d-01)
        parameter(d12h4l7 =  5.12309415013476d-01)
        parameter(d12h4l8 =  6.53804118732517d-01)
        parameter(d12h4l9 =  3.69560118619350d-01)
        parameter(d12h4l10 = -3.99395325110421d-02)
        parameter(d12h4l11 = -9.47945572764480d-02)
        parameter(d12h4l12 =  2.98655864523818d-03)
        parameter(d12h4l13 =  1.31793384194650d-02)
c
        parameter(d12h5l1 =  3.16485336386148d-02)
        parameter(d12h5l2 = -8.21447196104612d-02)
        parameter(d12h5l3 =  4.32755105849784d-02)
        parameter(d12h5l4 =  6.83065393224048d-02)
        parameter(d12h5l5 = -4.84004814948902d-03)
        parameter(d12h5l6 = -7.32083409088707d-02)
        parameter(d12h5l7 = -5.78691815689817d-02)
        parameter(d12h5l8 =  2.05585520869158d-02)
        parameter(d12h5l9 =  4.00194420635028d-01)
        parameter(d12h5l10 =  7.64498356903161d-01)
        parameter(d12h5l11 =  4.65380251184249d-01)
        parameter(d12h5l12 = -4.68642116294404d-02)
        parameter(d12h5l13 = -1.13713697920406d-01)
        parameter(d12h5l14 =  3.41329731746386d-03)
        parameter(d12h5l15 =  1.50624869010468d-02)
c
        parameter(d12h6l1 = -5.55080497108721d-03)
        parameter(d12h6l2 =  1.47108792588265d-02)
        parameter(d12h6l3 = -8.85911791140839d-03)
        parameter(d12h6l4 = -1.15892336474091d-02)
        parameter(d12h6l5 =  1.56983225887670d-03)
        parameter(d12h6l6 =  1.15711117119591d-02)
        parameter(d12h6l7 =  7.44552668496767d-03)
        parameter(d12h6l8 =  2.90249390308689d-02)
        parameter(d12h6l9 = -2.99994097480762d-02)
        parameter(d12h6l10 = -6.57634046784125d-02)
        parameter(d12h6l11 =  3.43825878587813d-01)
        parameter(d12h6l12 =  7.86713034635482d-01)
        parameter(d12h6l13 =  4.89584492012585d-01)
        parameter(d12h6l14 = -4.82546866951735d-02)
        parameter(d12h6l15 = -1.17799325128056d-01)
        parameter(d12h6l16 =  3.48847566515039d-03)
        parameter(d12h6l17 =  1.53942402679384d-02)
c
        parameter(d12g1l1 = -2.32170591862106d-02)
        parameter(d12g1l2 =  1.17136687081190d-01)
        parameter(d12g1l3 = -3.31407244748097d-01)
        parameter(d12g1l4 =  6.08338711236398d-01)
        parameter(d12g1l5 = -6.31390527628537d-01)
        parameter(d12g1l6 =  3.21364401968350d-01)
        parameter(d12g1l7 = -6.24831018002666d-02)
c
        parameter(d12g2l1 = -6.09230745955748d-02)
        parameter(d12g2l2 =  2.49618502573294d-01)
        parameter(d12g2l3 = -4.76418991643954d-01)
        parameter(d12g2l4 =  3.56579295733774d-01)
        parameter(d12g2l5 =  2.89304796889059d-01)
        parameter(d12g2l6 = -6.33606169200403d-01)
        parameter(d12g2l7 =  3.06976719889313d-01)
        parameter(d12g2l8 = -4.77281546675813d-03)
        parameter(d12g2l9 = -2.10618834993764d-02)
c
        parameter(d12g3l1 = -1.14630921543860d-01)
        parameter(d12g3l2 =  3.71547680503188d-01)
        parameter(d12g3l3 = -4.20756472301043d-01)
        parameter(d12g3l4 = -1.03448116555525d-01)
        parameter(d12g3l5 =  4.06710810431049d-01)
        parameter(d12g3l6 =  2.97735191981046d-01)
        parameter(d12g3l7 = -6.14859915730224d-01)
        parameter(d12g3l8 =  6.85311500641674d-02)
        parameter(d12g3l9 =  1.53810026669797d-01)
        parameter(d12g3l10 = -5.44887448013345d-03)
        parameter(d12g3l11 = -2.40452538554246d-02)
c
        parameter(d12g4l1 = -1.15086586235804d-01)
        parameter(d12g4l2 =  2.92452879358586d-01)
        parameter(d12g4l3 = -1.45671753991530d-01)
        parameter(d12g4l4 = -2.55255397973553d-01)
        parameter(d12g4l5 =  1.27975081910958d-03)
        parameter(d12g4l6 =  3.05194428781224d-01)
        parameter(d12g4l7 =  4.35229206763049d-01)
        parameter(d12g4l8 = -7.15563122827798d-01)
        parameter(d12g4l9 =  1.14554275482418d-01)
        parameter(d12g4l10 =  4.30036072497578d-02)
        parameter(d12g4l11 =  6.98774908402524d-02)
        parameter(d12g4l12 = -4.39592805179993d-03)
        parameter(d12g4l13 = -1.93987228593900d-02)
c
        parameter(d12g5l1 =  7.16575244096942d-02)
        parameter(d12g5l2 = -1.86977183443182d-01)
        parameter(d12g5l3 =  1.02112552332146d-01)
        parameter(d12g5l4 =  1.53385042989785d-01)
        parameter(d12g5l5 = -1.33046557036333d-02)
        parameter(d12g5l6 = -1.61626977522510d-01)
        parameter(d12g5l7 = -1.22226141624071d-01)
        parameter(d12g5l8 = -2.14298948728626d-01)
        parameter(d12g5l9 =  8.10373335613495d-01)
        parameter(d12g5l10 = -4.05884152396525d-01)
        parameter(d12g5l11 = -1.13504585558578d-01)
        parameter(d12g5l12 =  2.62163396939316d-02)
        parameter(d12g5l13 =  6.07171478855183d-02)
        parameter(d12g5l14 = -2.01560040142400d-03)
        parameter(d12g5l15 = -8.89461181387844d-03)
c
        parameter(d12g6l1 = -1.09269319891582d-02)
        parameter(d12g6l2 =  2.89588227291709d-02)
        parameter(d12g6l3 = -1.74394487657404d-02)
        parameter(d12g6l4 = -2.28137663872740d-02)
        parameter(d12g6l5 =  3.09026355933606d-03)
        parameter(d12g6l6 =  2.27780927944879d-02)
        parameter(d12g6l7 =  1.46567505314764d-02)
        parameter(d12g6l8 =  8.67342233103233d-02)
        parameter(d12g6l9 = -6.57618860628415d-02)
        parameter(d12g6l10 = -4.76210522536601d-01)
        parameter(d12g6l11 =  7.96862238947179d-01)
        parameter(d12g6l12 = -3.39909602687197d-01)
        parameter(d12g6l13 = -7.59730955669994d-02)
        parameter(d12g6l14 =  2.11811764748409d-02)
        parameter(d12g6l15 =  4.51380900148493d-02)
        parameter(d12g6l16 = -1.77212122147796d-03)
        parameter(d12g6l17 = -7.82016630927769d-03)
c
c
c        right edge coefficients for the wavelets on the interval
c
c
c
        parameter(d12h1r1 =  9.23667527528915d-01)
        parameter(d12h1r2 =  3.79334388520850d-01)
        parameter(d12h1r3 =  4.98645977823906d-02)
        parameter(d12h1r4 = -1.68894716290216d-02)
        parameter(d12h1r5 = -3.31180699766844d-03)
        parameter(d12h1r6 =  2.80462964022002d-03)
        parameter(d12h1r7 = -1.23755362396423d-02)
c
        parameter(d12h2r1 = -2.83195507191535d-01)
        parameter(d12h2r2 =  6.20923259598963d-01)
        parameter(d12h2r3 =  6.49443536714093d-01)
        parameter(d12h2r4 =  3.27901673168667d-01)
        parameter(d12h2r5 =  4.61346542276358d-02)
        parameter(d12h2r6 = -2.20110504191574d-02)
        parameter(d12h2r7 =  4.77299626789720d-02)
        parameter(d12h2r8 =  1.81108386509095d-03)
        parameter(d12h2r9 = -7.99210395621168d-03)
c
        parameter(d12h3r1 =  1.41588995000149d-01)
        parameter(d12h3r2 = -3.26411735373839d-01)
        parameter(d12h3r3 =  5.89871851747429d-02)
        parameter(d12h3r4 =  5.66257549421352d-01)
        parameter(d12h3r5 =  6.51143624454035d-01)
        parameter(d12h3r6 =  3.43470984054783d-01)
        parameter(d12h3r7 = -6.89756591213214d-02)
        parameter(d12h3r8 = -2.13985226417970d-02)
        parameter(d12h3r9 =  4.51454054268944d-02)
        parameter(d12h3r10 =  1.80701982947388d-03)
        parameter(d12h3r11 = -7.97416983634047d-03)
c
        parameter(d12h4r1 = -4.08092914580912d-02)
        parameter(d12h4r2 =  1.25197735301999d-01)
        parameter(d12h4r3 = -9.20894529766251d-02)
        parameter(d12h4r4 = -1.81544160431768d-01)
        parameter(d12h4r5 =  7.49289363228677d-02)
        parameter(d12h4r6 =  4.92189016559278d-01)
        parameter(d12h4r7 =  7.59888540158075d-01)
        parameter(d12h4r8 =  3.29853586406158d-01)
        parameter(d12h4r9 = -6.88520353723187d-02)
        parameter(d12h4r10 = -2.05622468853270d-02)
        parameter(d12h4r11 =  4.35703809518960d-02)
        parameter(d12h4r12 =  1.72945876176131d-03)
        parameter(d12h4r13 = -7.63190180112583d-03)
c
        parameter(d12h5r1 =  1.47333838526122d-02)
        parameter(d12h5r2 = -3.87992830632594d-02)
        parameter(d12h5r3 =  3.22268544105171d-02)
        parameter(d12h5r4 =  4.18135174610767d-02)
        parameter(d12h5r5 = -4.47016040254919d-02)
        parameter(d12h5r6 = -1.14359280112986d-01)
        parameter(d12h5r7 = -3.07926459677964d-02)
        parameter(d12h5r8 =  4.95330334080175d-01)
        parameter(d12h5r9 =  7.83106303708429d-01)
        parameter(d12h5r10 =  3.36526176854865d-01)
        parameter(d12h5r11 = -7.16155454146472d-02)
        parameter(d12h5r12 = -2.09664597972113d-02)
        parameter(d12h5r13 =  4.44599515226543d-02)
        parameter(d12h5r14 =  1.76224419266212d-03)
        parameter(d12h5r15 = -7.77658012169345d-03)
c
        parameter(d12h6r1 = -2.27977846204926d-03)
        parameter(d12h6r2 =  5.93426712836801d-03)
        parameter(d12h6r3 = -4.65392062045452d-03)
        parameter(d12h6r4 = -4.69732619853451d-03)
        parameter(d12h6r5 =  7.90549825302031d-03)
        parameter(d12h6r6 =  1.47842676156550d-02)
        parameter(d12h6r7 =  6.35660984178306d-04)
        parameter(d12h6r8 = -1.18413762192955d-01)
        parameter(d12h6r9 = -4.71952433748829d-02)
        parameter(d12h6r10 =  4.91405649375928d-01)
        parameter(d12h6r11 =  7.87382709556111d-01)
        parameter(d12h6r12 =  3.37877716874012d-01)
        parameter(d12h6r13 = -7.25679362426239d-02)
        parameter(d12h6r14 = -2.10569052473585d-02)
        parameter(d12h6r15 =  4.47137346291320d-02)
        parameter(d12h6r16 =  1.76757325116928d-03)
        parameter(d12h6r17 = -7.80009664149626d-03)
c
        parameter(d12g1r1 =  5.56597135210399d-02)
        parameter(d12g1r2 = -2.24249808504556d-01)
        parameter(d12g1r3 =  5.13542423752636d-01)
        parameter(d12g1r4 = -6.32191648038890d-01)
        parameter(d12g1r5 =  4.90183861475066d-01)
        parameter(d12g1r6 = -2.04232228314738d-01)
        parameter(d12g1r7 =  3.50841516864278d-02)
c
        parameter(d12g2r1 = -1.05715927784490d-01)
        parameter(d12g2r2 =  3.39436457392220d-01)
        parameter(d12g2r3 = -4.78785450382155d-01)
        parameter(d12g2r4 =  7.18764547510599d-02)
        parameter(d12g2r5 =  5.04278061448007d-01)
        parameter(d12g2r6 = -5.79696670464663d-01)
        parameter(d12g2r7 =  2.20535578330557d-01)
        parameter(d12g2r8 =  3.31337148263599d-03)
        parameter(d12g2r9 = -1.46215257311921d-02)
c
        parameter(d12g3r1 =  1.74329041850787d-01)
        parameter(d12g3r2 = -4.25553144657815d-01)
        parameter(d12g3r3 =  2.63804709985222d-01)
        parameter(d12g3r4 =  3.61643820384694d-01)
        parameter(d12g3r5 = -2.64467105509555d-01)
        parameter(d12g3r6 = -4.90854706267411d-01)
        parameter(d12g3r7 =  4.96210757800414d-01)
        parameter(d12g3r8 =  6.82842769425827d-02)
        parameter(d12g3r9 = -1.60689947480903d-01)
        parameter(d12g3r10 = -5.15667378981358d-03)
        parameter(d12g3r11 =  2.27558059518092d-02)
c
        parameter(d12g4r1 =  3.24470951782164d-03)
        parameter(d12g4r2 =  5.11125346019136d-04)
        parameter(d12g4r3 =  4.38149222239852d-03)
        parameter(d12g4r4 = -3.03584706528650d-02)
        parameter(d12g4r5 = -3.79693344330058d-02)
        parameter(d12g4r6 =  7.32218690931611d-02)
        parameter(d12g4r7 =  3.43681870756667d-01)
        parameter(d12g4r8 = -7.86742846062982d-01)
        parameter(d12g4r9 =  4.88619040150674d-01)
        parameter(d12g4r10 =  4.82381115811522d-02)
        parameter(d12g4r11 = -1.17637716587594d-01)
        parameter(d12g4r12 = -3.49171926751491d-03)
        parameter(d12g4r13 =  1.54085539106081d-02)
c
        parameter(d12g5r1 =  5.21103425503510d-03)
        parameter(d12g5r2 = -1.37905434569141d-02)
        parameter(d12g5r3 =  1.17228748724916d-02)
        parameter(d12g5r4 =  1.65183964462511d-02)
        parameter(d12g5r5 = -1.48460330074742d-02)
        parameter(d12g5r6 = -4.32876920273277d-02)
        parameter(d12g5r7 = -1.49192354728936d-02)
        parameter(d12g5r8 =  7.35272868470368d-02)
        parameter(d12g5r9 =  3.35344107657852d-01)
        parameter(d12g5r10 = -7.88017770953177d-01)
        parameter(d12g5r11 =  4.91636285183262d-01)
        parameter(d12g5r12 =  4.83280932891972d-02)
        parameter(d12g5r13 = -1.18025958966888d-01)
        parameter(d12g5r14 = -3.49204326328664d-03)
        parameter(d12g5r15 =  1.54099836665344d-02)
c
        parameter(d12g6r1 = -1.15437604256422d-03)
        parameter(d12g6r2 =  3.00484275871557d-03)
        parameter(d12g6r3 = -2.35653356573041d-03)
        parameter(d12g6r4 = -2.37851217474146d-03)
        parameter(d12g6r5 =  4.00298447403386d-03)
        parameter(d12g6r6 =  7.48608017246960d-03)
        parameter(d12g6r7 =  3.21869788465603d-04)
        parameter(d12g6r8 = -4.49338348826036d-02)
        parameter(d12g6r9 = -2.04926062854434d-02)
        parameter(d12g6r10 =  7.27936391585918d-02)
        parameter(d12g6r11 =  3.37760492983485d-01)
        parameter(d12g6r12 = -7.87664903147568d-01)
        parameter(d12g6r13 =  4.91084187805316d-01)
        parameter(d12g6r14 =  4.83133225841471d-02)
        parameter(d12g6r15 = -1.17995169148154d-01)
        parameter(d12g6r16 = -3.49078227503662d-03)
        parameter(d12g6r17 =  1.54044190710035d-02)
c
c
c        left edge coefficients for the preconditioner
c
c
c
        parameter(d12pc1l1 =  2.69125823813972d+00)
        parameter(d12pc1l2 = -4.69139083476636d+00)
        parameter(d12pc2l2 =  2.27066616576605d+00)
        parameter(d12pc1l3 =  5.65087880986530d+00)
        parameter(d12pc2l3 = -2.34621803490825d+00)
        parameter(d12pc3l3 =  1.60292673096433d+00)
        parameter(d12pc1l4 = -4.00361022496952d+00)
        parameter(d12pc2l4 =  1.68696449245117d+00)
        parameter(d12pc3l4 = -7.48546437672687d-01)
        parameter(d12pc4l4 =  1.16880748007002d+00)
        parameter(d12pc1l5 =  1.54810766871523d+00)
        parameter(d12pc2l5 = -6.55773003218562d-01)
        parameter(d12pc3l5 =  3.03044772965090d-01)
        parameter(d12pc4l5 = -1.28788492546969d-01)
        parameter(d12pc5l5 =  1.02268034676095d+00)
        parameter(d12pc1l6 = -2.52340379611320d-01)
        parameter(d12pc2l6 =  1.06944336176567d-01)
        parameter(d12pc3l6 = -5.08302032095054d-02)
        parameter(d12pc4l6 =  2.18056797614183d-02)
        parameter(d12pc5l6 = -7.65254204100803d-03)
        parameter(d12pc6l6 =  1.00064108776624d+00)
c
        parameter(d12ipc1l1 =  3.71573409726460d-01)
        parameter(d12ipc1l2 =  7.67702498550905d-01)
        parameter(d12ipc2l2 =  4.40399392511595d-01)
        parameter(d12ipc1l3 = -1.86233627540046d-01)
        parameter(d12ipc2l3 =  6.44616486401546d-01)
        parameter(d12ipc3l3 =  6.23858833147287d-01)
        parameter(d12ipc1l4 =  4.54683334273958d-02)
        parameter(d12ipc2l4 = -2.22802101752090d-01)
        parameter(d12ipc3l4 =  3.99541682549008d-01)
        parameter(d12ipc4l4 =  8.55572895495239d-01)
        parameter(d12ipc1l5 = -9.29336969083579d-03)
        parameter(d12ipc2l5 =  6.33238224006255d-02)
        parameter(d12ipc3l5 = -1.34549165713442d-01)
        parameter(d12ipc4l5 =  1.07744266156934d-01)
        parameter(d12ipc5l5 =  9.77822643377489d-01)
        parameter(d12ipc1l6 =  1.13193655198438d-03)
        parameter(d12ipc2l6 = -8.98353496027236d-03)
        parameter(d12ipc3l6 =  2.19548751328578d-02)
        parameter(d12ipc4l6 = -1.78204066006140d-02)
        parameter(d12ipc5l6 =  7.47803481046330d-03)
        parameter(d12ipc6l6 =  9.99359322963966d-01)
c
c
c        right edge coefficients for the preconditioner
c
c
c
        parameter(d12pc1r1 =  6.30332954789026d-01)
        parameter(d12pc1r2 =  1.50459696638497d-01)
        parameter(d12pc2r2 =  9.76051909205893d-01)
        parameter(d12pc1r3 =  1.46150031713424d-01)
        parameter(d12pc2r3 = -7.02619642126096d-02)
        parameter(d12pc3r3 =  9.78247075887753d-01)
        parameter(d12pc1r4 = -1.80473694617946d-01)
        parameter(d12pc2r4 =  1.47937536800329d-01)
        parameter(d12pc3r4 = -2.39282102905875d-02)
        parameter(d12pc4r4 =  1.02211853982210d+00)
        parameter(d12pc1r5 =  7.50231007888198d-02)
        parameter(d12pc2r5 = -6.71131275252377d-02)
        parameter(d12pc3r5 =  9.36848497496388d-03)
        parameter(d12pc4r5 = -1.66413898429385d-02)
        parameter(d12pc5r5 =  1.00310267533560d+00)
        parameter(d12pc1r6 = -1.24803281019971d-02)
        parameter(d12pc2r6 =  1.19527720534302d-02)
        parameter(d12pc3r6 = -1.42903840021121d-03)
        parameter(d12pc4r6 =  2.63346846307428d-03)
        parameter(d12pc5r6 = -9.85097652640761d-04)
        parameter(d12pc6r6 =  1.00007841999455d+00)
c
        parameter(d12ipc1r1 =  1.58646314206228d+00)
        parameter(d12ipc1r2 = -2.44555397957319d-01)
        parameter(d12ipc2r2 =  1.02453567332663d+00)
        parameter(d12ipc1r3 = -2.54582494834295d-01)
        parameter(d12ipc2r3 =  7.35866128181275d-02)
        parameter(d12ipc3r3 =  1.02223663596695d+00)
        parameter(d12ipc1r4 =  3.09555175851304d-01)
        parameter(d12ipc2r4 = -1.46564690976056d-01)
        parameter(d12ipc3r4 =  2.39309749693193d-02)
        parameter(d12ipc4r4 =  9.78360103099246d-01)
        parameter(d12ipc1r5 = -1.27502183298440d-01)
        parameter(d12ipc2r5 =  6.54283551181634d-02)
        parameter(d12ipc3r5 = -9.15017386236866d-03)
        parameter(d12ipc4r5 =  1.62309126301601d-02)
        parameter(d12ipc5r5 =  9.96906921482821d-01)
        parameter(d12ipc1r6 =  2.14164019309600d-02)
        parameter(d12ipc2r6 = -1.16895397642489d-02)
        parameter(d12ipc3r6 =  1.38867122464423d-03)
        parameter(d12ipc4r6 = -2.56029066512977d-03)
        parameter(d12ipc5r6 =  9.81973661884842d-04)
        parameter(d12ipc6r6 =  9.99921586154663d-01)
c