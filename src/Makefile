MARSOBJECTS = dcalcvar.o dmarss.o dorthreg.o dqrreg.o S_compat.o
BRUTOOBJECTS = mspline.o bruto.o sortdi.o S_compat.o


all: MARS.so

MARS.so: $(MARSOBJECTS)
	$(LD) $(SHLIBLDFLAGS) -o MARS.so $(MARSOBJECTS)

BRUTO.so: $(BRUTOOBJECTS)
	$(LD) $(SHLIBLDFLAGS) -o BRUTO.so $(BRUTOOBJECTS)

clean:
	rm -f *.o
	rm -f *.so

realclean: clean


