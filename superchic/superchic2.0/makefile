LIBFLAG = LHAPDF
LHAPDFLIB = /home/lang/lhapdf614/lib

#####################

HOME = $(PWD)
MAIN = $(PWD)/src/main
SOURCEDIR = $(PWD)/src
VPATH = $(DIRS)
INCPATH = $(SOURCEDIR)/inc

FC = gfortran
FFLAGS 	= -fno-automatic -fno-f2c -O2 -g -I$(INCPATH)	

OUTPUT_OPTION = -o $(HOME)/obj/$@

DIRS	 =	$(SOURCEDIR)/int:\
		$(SOURCEDIR)/main:\
		$(SOURCEDIR)/mes:\
		$(SOURCEDIR)/PDFs:\
		$(SOURCEDIR)/phase:\
		$(SOURCEDIR)/sPDFs:\
		$(SOURCEDIR)/subamps:\
		$(SOURCEDIR)/surv:\
		$(SOURCEDIR)/EW:\
		$(SOURCEDIR)/unw:\
		$(SOURCEDIR)/user:\
		$(SOURCEDIR)/int:\
		$(SOURCEDIR)/var:\
		$(SOURCEDIR)/init:\

#############

PDFINPUT     = 	USER

#############

Mesf = \
calcmes.o \
mesint.o \
wfinit.o \
wfoctet.o \
wfsinglet.o \

PDFsfLHA = \
alphas.o \
inpdf.o \

PDFsfUSER = \
alphasuser.o \
inpdfuser.o \


Intf = \
ran.o  \
vegas.o \

Phasef = \
2body.o \
2bodyw.o \
2jetps.o \
2jetpsm.o \
3jetps.o \
boost.o \
chic0decay3.o \
chic1decay3.o \
chic1decay2s.o \
chic1decay2f.o \
chic2decay3.o \
chic2decay2s.o \
chic2decay2f.o \
jpsidecayphot.o \
genpol1.o \
genpol2.o \
rambo.o \
3body.o \
3bodyinit.o \
2bodyinit.o \
wwcorr.o \
jpsidecay.o \

Subampsf = \
chi0.o \
chi1.o \
chi2.o \
etaq.o \
higgs.o \
pipi.o \
qqjets.o \
diphoton.o \
etaeta.o \
etapetap.o \
etaetap.o \
eta.o \
gggjets.o \
pipixy.o \
rhorho.o \
djpsi.o \
djpsip.o \
djpsipp.o \
ggjets.o \
qqgjets.o \
rhorhoxy.o \
ww.o \
wwpol.o \
ll.o \
llpol.o \
mhv.o \
lightlight.o \
lightlightpol.o \

Survf = \
initparsr.o \
formfac.o \
formfacphot.o \
formfacgam.o \
seik.o \
seikphot.o\
seikgam.o \
screeningint.o \
readscreen.o \
formfacgamel.o \

Userf = \
cuts.o \
histo.o \

Mainf = \
bare.o \
header.o \
main.o \
process.o \
superchic.o \
wtgen.o \

sPDFsf = \
calchg.o \
hpdfint.o \
calcsud.o \
sudint.o \
sPDF.o \

Unwf = \
unweight.o \
unwprint.o \

Varf = \
mu.o \
nf.o \
string.o \
varfuncs.o \

InitfLHA = \
alphas.o \
init.o \
initsud.o \
nf.o \
string.o \
hg.o \
inithg.o \
initpars.o \
calcop.o \
calcscreen.o \
opacityint.o \
screeningint.o \
screening.o \
opacity.o \
PDF.o \
PDFlha.o \
Sudakov.o \
inpdf.o \

InitfUSER = \
alphasuser.o \
init.o \
initsud.o \
nf.o \
string.o \
hg.o \
inithg.o \
initpars.o \
calcop.o \
calcscreen.o \
opacityint.o \
screeningint.o \
screening.o \
opacity.o \
PDF.o \
PDFuser.o \
Sudakov.o \
inpdfuser.o \


iCODELHA = $(InitfLHA) \

sCODELHA = $(Mainf) $(Mesf) $(EW) $(PDFsfLHA) $(Intf) \
$(Phasef) $(Subampsf) $(Survf) $(Userf) \
$(sPDFsf) $(Unwf) $(Varf) \

iCODEUSER = $(InitfUSER) \

sCODEUSER = $(Mainf) $(Mesf) $(EW) $(PDFsfUSER) $(Intf) \
$(Phasef) $(Subampsf) $(Survf) $(Userf) \
$(sPDFsf) $(Unwf) $(Varf) \

clean: 	
	-rm -f $(HOME)/obj/*.o
	make all

all: init superchic

ifeq ($(PDFINPUT),LHAPDF)
init:	$(iCODELHA)
	$(FC) $(FFLAGS) -Wl,-R$(LHAPDFLIB) -L$(LHAPDFLIB) -l$(LIBFLAG) -o $@ \
	$(patsubst %,obj/%,$(iCODELHA))  -l$(LIBFLAG) 
	mv init bin/init
	@echo '    ----> Init compiled OK <----'

superchic:	$(sCODELHA)
	$(FC) $(FFLAGS) -Wl,-R$(LHAPDFLIB) -L$(LHAPDFLIB) -o $@ \
	$(patsubst %,obj/%,$(sCODELHA)) -l$(LIBFLAG)
	mv superchic bin/superchic
	@echo '    ----> Superchicv2.0 compiled OK <----'
endif

ifeq ($(PDFINPUT),USER)
init:	$(iCODEUSER)
	$(FC) $(FFLAGS)  -o $@ \
	$(patsubst %,obj/%,$(iCODEUSER))
	mv init bin/init
	@echo '    ----> Init compiled OK <----'

superchic:	$(sCODEUSER)
	$(FC) $(FFLAGS)  -o $@ \
	$(patsubst %,obj/%,$(sCODEUSER)) 
	mv superchic bin/superchic
	@echo '    ----> Superchicv2.0 compiled OK <----'
endif
