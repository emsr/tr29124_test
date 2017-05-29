c	$Log:	D08qmf.f,v $
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
c        Daubechies' Interval QMF's with lf = 8                             c
c                                                                           c
cvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvc
c
         parameter(d8h1 = 3.22231006040514678716159d-02)
         parameter(d8h2 =-1.26039672620313037539160d-02)
         parameter(d8h3 =-9.92195435766335325852080d-02)
         parameter(d8h4 = 2.97857795605306051402901d-01)
         parameter(d8h5 = 8.03738751805132080878805d-01)
         parameter(d8h6 = 4.97618667632774989979605d-01)
         parameter(d8h7 =-2.96355276460024917643691d-02)
         parameter(d8h8 =-7.57657147895022132277461d-02)
c
         parameter(d8g1 = d8h8)
         parameter(d8g2 =-d8h7)
         parameter(d8g3 = d8h6)
         parameter(d8g4 =-d8h5)
         parameter(d8g5 = d8h4)
         parameter(d8g6 =-d8h3)
         parameter(d8g7 = d8h2)
         parameter(d8g8 =-d8h1)
c
c        left edge coefficients for the wavelets on the interval
c
        parameter(d8h1l1 = .90975392578299541422d0)
        parameter(d8h1l2 = .40416588940201723557d0)
        parameter(d8h1l3 = .089040318656793899574d0)
        parameter(d8h1l4 = -.011984192014966211145d0)
        parameter(d8h1l5 = -.030429084139182777526d0)
c
        parameter(d8g1l1 = .075739707618965346219d0)
        parameter(d8g1l2 = -.32543917184495103305d0)
        parameter(d8g1l3 = .68434908047295674723d0)
        parameter(d8g1l4 = -.62004421070554711214d0)
        parameter(d8g1l5 = .18858513977781955539d0)
c
        parameter(d8h2l1 = -.27285140766949225591d0)
        parameter(d8h2l2 = .50908152323149266343d0)
        parameter(d8h2l3 = .62364244337236631118d0)
        parameter(d8h2l4 = .46284008632544893219d0)
        parameter(d8h2l5 = .24674764172786638137d0)
        parameter(d8h2l6 = -.017669533291129221108d0)
        parameter(d8h2l7 = -.045173645490327315413d0)
c
        parameter(d8g2l1 = .16659595920609478902d0)
        parameter(d8g2l2 = -.48478430891405802583d0)
        parameter(d8g2l3 = .35646354251378829101d0)
        parameter(d8g2l4 = .48398961557429770139d0)
        parameter(d8g2l5 = -.60575436513594114671d0)
        parameter(d8g2l6 = .034518332894524950210d0)
        parameter(d8g2l7 = .088249016394632876783d0)
c
        parameter(d8h3l1 = .12611792859809708119d0)
        parameter(d8h3l2 = -.23085572680672168167d0)
        parameter(d8h3l3 = -.052799235246625577463d0)
        parameter(d8h3l4 = .21926517128941860346d0)
        parameter(d8h3l5 = .46348072109960973065d0)
        parameter(d8h3l6 = .70011971404312468766d0)
        parameter(d8h3l7 = .41203257897604627778d0)
        parameter(d8h3l8 = -.026222762497916067111d0)
        parameter(d8h3l9 = -.067040694133818099667d0)
c
        parameter(d8g3l1 = .20825353261585996816d0)
        parameter(d8g3l2 = -.40182279321111657707d0)
        parameter(d8g3l3 = -.068721487621235874441d0)
        parameter(d8g3l4 = .33021351133522513671d0)
        parameter(d8g3l5 = .55802129618734586368d0)
        parameter(d8g3l6 = -.59949741336936802557d0)
        parameter(d8g3l7 = -.069091991978684908007d0)
        parameter(d8g3l8 = .027853569971929078366d0)
        parameter(d8g3l9 = .071209990372730355011d0)
c
        parameter(d8h4l1 = -.029079804272549243429d0)
        parameter(d8h4l2 = .059928072294318161380d0)
        parameter(d8h4l3 = .0061764277780953420663d0)
        parameter(d8h4l4 = -.040210999043085892199d0)
        parameter(d8h4l5 = -.039525870128122906889d0)
        parameter(d8h4l6 = -.052599062573631102042d0)
        parameter(d8h4l7 = .32894944796959845562d0)
        parameter(d8h4l8 = .79663789672531588080d0)
        parameter(d8h4l9 = .49011303364021253823d0)
        parameter(d8h4l10 = -.029432877677387671779d0)
        parameter(d8h4l11 = -.075247623129128381429d0)
c
        parameter(d8g4l1 = -.065485007015422398501d0)
        parameter(d8g4l2 = .13495242945354731343d0)
        parameter(d8g4l3 = .013908739295079390932d0)
        parameter(d8g4l4 = -.090551419457775651279d0)
        parameter(d8g4l5 = -.089008573041674695765d0)
        parameter(d8g4l6 = -.37334444756195190139d0)
        parameter(d8g4l7 = .84046537082509037751d0)
        parameter(d8g4l8 = -.31568493615246357870d0)
        parameter(d8g4l9 = -.12029765089741391514d0)
        parameter(d8g4l10 = .013070202799775884201d0)
        parameter(d8g4l11 = .033415070904004976114d0)
c
c
c
c        right edge coefficients for the wavelets on the interval
c
        parameter(d8h1r1 = .91547051883810854018d0)
        parameter(d8h1r2 = .39191428098156657169d0)
        parameter(d8h1r3 = .059477711238407501787d0)
        parameter(d8h1r4 = -.025191808506643957308d0)
        parameter(d8h1r5 = .064379345686262128356d0)
c
        parameter(d8g1r1 = .19827799059709502246d0)
        parameter(d8g1r2 = -.60406778678838103018d0)
        parameter(d8g1r3 = .64952976030516477709d0)
        parameter(d8g1r4 = -.40503095407943504406d0)
        parameter(d8g1r5 = .099241947405233368263d0)
c
        parameter(d8h2r1 = -.21916264686360397441d0)
        parameter(d8h2r2 = .44880017813063045390d0)
        parameter(d8h2r3 = .75400050839596823140d0)
        parameter(d8h2r4 = .39377581566148591390d0)
        parameter(d8h2r5 = -.15813389438931389654d0)
        parameter(d8h2r6 = -.016142011904287920593d0)
        parameter(d8h2r7 = .041268408805739583103d0)
c
        parameter(d8g2r1 = -.27262735058992725801d0)
        parameter(d8g2r2 = .50928674835762860473d0)
        parameter(d8g2r3 = .068118566264738166854d0)
        parameter(d8g2r4 = -.67353457518648896953d0)
        parameter(d8g2r5 = .44993408265114559900d0)
        parameter(d8g2r6 = .027190655403678115961d0)
        parameter(d8g2r7 = -.069515193617030161473d0)
c
        parameter(d8h3r1 = .012900782893843634354d0)
        parameter(d8h3r2 = -.13907160056079416442d0)
        parameter(d8h3r3 = .029213679496076078323d0)
        parameter(d8h3r4 = .46061685367679247208d0)
        parameter(d8h3r5 = .81641197419298729535d0)
        parameter(d8h3r6 = .29864733458521474811d0)
        parameter(d8h3r7 = -.10276635357182394708d0)
        parameter(d8h3r8 = -.012574882110017058721d0)
        parameter(d8h3r9 = .032148741970777129858d0)
c
        parameter(d8g3r1 = -.0045819593397576190533d0)
        parameter(d8g3r2 = -.030620323648208387712d0)
        parameter(d8g3r3 = -.013887371452971361363d0)
        parameter(d8g3r4 = .095048349988442514248d0)
        parameter(d8g3r5 = .30158150124129352334d0)
        parameter(d8g3r6 = -.80336363165658830946d0)
        parameter(d8g3r7 = .49684270650616092637d0)
        parameter(d8g3r8 = .029632043697507898216d0)
        parameter(d8g3r9 = -.075756807782644236930d0)
c
        parameter(d8h4r1 = -.0067756036517274764721d0)
        parameter(d8h4r2 = .019132441289899007282d0)
        parameter(d8h4r3 = -.017709184245907807129d0)
        parameter(d8h4r4 = -.067659161740038464179d0)
        parameter(d8h4r5 = -.030235884812807082946d0)
        parameter(d8h4r6 = .49779208211653976270d0)
        parameter(d8h4r7 = .80394959961412010659d0)
        parameter(d8h4r8 = .29771110110016361257d0)
        parameter(d8h4r9 = -.099108040550099571666d0)
        parameter(d8h4r10 = -.012598951895850180172d0)
        parameter(d8h4r11 = .032210278399291594766d0)
c
        parameter(d8g4r1 = -.0028803051241456110658d0)
        parameter(d8g4r2 = .0081331895307455420533d0)
        parameter(d8g4r3 = -.0075281637990915038805d0)
        parameter(d8g4r4 = -.028761869830674646117d0)
        parameter(d8g4r5 = -.012853256836710187824d0)
        parameter(d8g4r6 = .099201914392120836400d0)
        parameter(d8g4r7 = .29778998412792446390d0)
        parameter(d8g4r8 = -.80379368413398200528d0)
        parameter(d8g4r9 = .49764705235425275855d0)
        parameter(d8g4r10 = .029637660176292392689d0)
        parameter(d8g4r11 = -.075771166782247360195d0)
c
c        left edge coefficients for the preconditioner
c
c
	 parameter(d8pc1l1 = 2.4899111140824833432d0)
	 parameter(d8pc1l2 =-2.7529885357135795429d0)
	 parameter(d8pc2l2 = 1.6772105498044283777d0)
	 parameter(d8pc1l3 = 1.6878414467970039542d0)
	 parameter(d8pc2l3 =-.70753754370482426899d0)
	 parameter(d8pc3l3 = 1.1301451419680759567d0)
	 parameter(d8pc1l4 =-.40222211735395905367d0)
	 parameter(d8pc2l4 = .17635442878554937940d0)
	 parameter(d8pc3l4 =-.061621215501286551203d0)
	 parameter(d8pc4l4 = 1.0068851564850727934d0)
c
	 parameter(d8ipc1l1 =  .40162076242167132109d0)
	 parameter(d8ipc1l2 =  .65922394465043988597d0)
	 parameter(d8ipc2l2 = .59622806457818029458d0)
	 parameter(d8ipc1l3 =-.18709674563738923634d0)
	 parameter(d8ipc2l3 = .37327394918930106812d0)
	 parameter(d8ipc3l3 = .88484209935952433399d0)
	 parameter(d8ipc1l4 = .033523746113502769821d0)
	 parameter(d8ipc2l4 =-.081584145680874628305d0)
	 parameter(d8ipc3l4 = .054152199322895070106d0)
	 parameter(d8ipc4l4 = .99316192473319585900d0)
c
c
c        right edge coefficients for the preconditioner
c
c
	 parameter(d8pc1r1 = .50051923113793834074d0)
	 parameter(d8pc1r2 = .37673863695833614763d0)
	 parameter(d8pc2r2 = .78081761658738585219d0)
	 parameter(d8pc1r3 =-.00093100685477696553087d0)
	 parameter(d8pc2r3 = .091704627903655303768d0)
	 parameter(d8pc3r3 = 1.0023129562376633347d0)
	 parameter(d8pc1r4 =-.0073733049057012719128d0)
	 parameter(d8pc2r4 =-.018445047273518395674d0)
	 parameter(d8pc3r4 =-.0022411542961159204922d0)
	 parameter(d8pc4r4 = 1.0003980780482839635d0)
c
	 parameter(d8ipc1r1 = 1.9979252300186034387d0)
	 parameter(d8ipc1r2 =-.96398392135616003463d0)
	 parameter(d8ipc2r2 = 1.2807088092742643304d0)
	 parameter(d8ipc1r3 = .090053578910487419586d0)
	 parameter(d8ipc2r3 =-.11717590207382437902d0)
	 parameter(d8ipc3r3 = .99769238118367204024d0)
	 parameter(d8ipc1r4 =-.0028464600220995317358d0)
	 parameter(d8ipc2r4 = .023350829801588023229d0)
	 parameter(d8ipc3r4 = .0022350928249024384967d0)
	 parameter(d8ipc4r4 = .99960208035479177484d0)
c
