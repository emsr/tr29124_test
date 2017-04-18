.SUFFIXES : .lib .ind

.PHONY : lib colors clean

lib : CFSF.mpl
	@maple -s -q $<

colors : colors.mpl
	@maple -s -q $<

clean :
	rm CFSF.lib CFSF.ind

