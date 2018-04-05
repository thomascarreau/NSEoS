CC = gcc
RM = rm
LIBDIR = ../nseos
BUILDDIR = ../nseos/build
OBJS = f_icrust4d.o \
	   $(BUILDDIR)/observables.o \
	   $(BUILDDIR)/nuclear_en.o \
	   $(BUILDDIR)/coulomb_en.o \
	   $(BUILDDIR)/nuclear_surface_en.o \
	   $(BUILDDIR)/nuclear_matter.o \
	   $(BUILDDIR)/empirical.o 
APP = icrust4d

$(APP) : $(OBJS) icrust4d.c
	$(CC) -o $(APP) icrust4d.c $(OBJS) -Wall -lgsl -lgslcblas -lm

$(BUILDDIR)/empirical.o : $(LIBDIR)/empirical.h \
	$(LIBDIR)/empirical.c 
	$(CC) -c $(LIBDIR)/empirical.c -o $(BUILDDIR)/empirical.o

$(BUILDDIR)/nuclear_matter.o : $(LIBDIR)/nuclear_matter.h \
	$(LIBDIR)/nuclear_matter.c \
	$(LIBDIR)/mathconst.h \
	$(LIBDIR)/phyconst.h
	$(CC) -c $(LIBDIR)/nuclear_matter.c -o $(BUILDDIR)/nuclear_matter.o

$(BUILDDIR)/nuclear_surface_en.o : $(LIBDIR)/nuclear_surface_en.h \
	$(LIBDIR)/nuclear_surface_en.c
	$(CC) -c $(LIBDIR)/nuclear_surface_en.c -o $(BUILDDIR)/nuclear_surface_en.o

$(BUILDDIR)/coulomb_en.o : $(LIBDIR)/coulomb_en.h \
	$(LIBDIR)/coulomb_en.c \
	$(LIBDIR)/empirical.h \
	$(LIBDIR)/mathconst.h \
	$(LIBDIR)/phyconst.h
	$(CC) -c $(LIBDIR)/coulomb_en.c -o $(BUILDDIR)/coulomb_en.o

$(BUILDDIR)/nuclear_en.o : $(LIBDIR)/nuclear_en.h \
	$(LIBDIR)/nuclear_en.c
	$(CC) -c $(LIBDIR)/nuclear_en.c -o $(BUILDDIR)/nuclear_en.o

$(BUILDDIR)/observables.o : $(LIBDIR)/observables.h \
	$(LIBDIR)/observables.c
	$(CC) -c $(LIBDIR)/observables.c -o $(BUILDDIR)/observables.o

f_icrust4d.o : f_icrust4d.h \
	f_icrust4d.c
	$(CC) -c f_icrust4d.c -o f_icrust4d.o

clean:
	$(RM) $(OBJS) $(APP)
