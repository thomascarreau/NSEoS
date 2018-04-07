CC = gcc
RM = rm -r
LIBDIR = ../nseos
MKDIR = mkdir -p
BUILDDIR = ../nseos/build
OBJS = f_core.o \
	   f_icrust4d.o \
	   $(BUILDDIR)/observables.o \
	   $(BUILDDIR)/nuclear_en.o \
	   $(BUILDDIR)/coulomb_en.o \
	   $(BUILDDIR)/nuclear_surface_en.o \
	   $(BUILDDIR)/nuclear_matter.o \
	   $(BUILDDIR)/empirical.o 
APP = nseos

$(APP) : $(OBJS) nseos.c
	$(CC) -o $(APP) nseos.c $(OBJS) -lgsl -lgslcblas -lm

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)

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

$(BUILDDIR)/observables.o : $(BUILDDIR) \
	$(LIBDIR)/observables.h \
	$(LIBDIR)/observables.c
	$(CC) -c $(LIBDIR)/observables.c -o $(BUILDDIR)/observables.o

f_icrust4d.o : f_icrust4d.h \
	f_icrust4d.c \
	choices.h
	$(CC) -c f_icrust4d.c -o f_icrust4d.o

f_core.o : f_core.h \
	f_core.c \
	choices.h
	$(CC) -c f_core.c -o f_core.o

clean:
	$(RM) $(OBJS) $(APP)
