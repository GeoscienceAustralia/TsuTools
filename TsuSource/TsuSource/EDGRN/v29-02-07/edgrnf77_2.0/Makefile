## Adding fortran compiler information
FC = gfortran
FFLAGS =	-O3


FORFILES =  axb.f       edgmain.f    bessj.f     edgsh.f      bessj0.f    \
            edgkern.f   edgmoment.f  bessj1.f    edgmatinv.f  edgsublay.f \
	    gemp.f      edgmatrix.f  getdata.f   edghask.f    edgbstab.f  \
            memcpy.f    edgwvint.f   edglayer.f  edgpsv.f

OBJECTS =	$(FORFILES:%.f=%.o)

PROGRAM =	edgrn2.0

$(PROGRAM): 	$(OBJECTS)
		$(FC) $(FFLAGS) $(OBJECTS) -o $@

bessj.o:bessj.f
axb.o:axb.f
edglayer.o:edgglobal.h edglayer.f
edgmatrix.o:edgglobal.h edgmatrix.f
bessj0.o:bessj0.f
bessj1.o:bessj1.f
edgwvint.o:edgglobal.h edgwvint.f
edgpsv.o:edgglobal.h edgpsv.f
edgmatinv.o:edgglobal.h edgmatinv.f
edgsublay.o:edgglobal.h edgsublay.f
edgmoment.o:edgglobal.h edgmoment.f
edgmain.o:edgglobal.h edgmain.f
memcpy.o:memcpy.f
gemp.o:gemp.f
edgbstab.o:edgglobal.h edgbstab.f
edgkern.o:edgglobal.h edgkern.f
edgsh.o:edgglobal.h edgsh.f
edghask.o:edgglobal.h edghask.f
getdata.o:getdata.f







