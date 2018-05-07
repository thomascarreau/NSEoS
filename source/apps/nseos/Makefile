CC = gcc
RM = rm -r
LIBDIR = ../../nseos
MKDIR = mkdir -p
BUILDDIR = ../../nseos/build
OBJS = $(BUILDDIR)/core.o \
	   $(BUILDDIR)/crust.o \
	   $(BUILDDIR)/nuclear_en.o \
	   $(BUILDDIR)/coulomb.o \
	   $(BUILDDIR)/nuclear_surface_en.o \
	   $(BUILDDIR)/nuclear_matter.o \
	   $(BUILDDIR)/empirical.o 
APP = nseos

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	nseos.c
	$(CC) -o $(APP) nseos.c $(OBJS) -lgsl -lgslcblas -lm

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)

$(BUILDDIR)/empirical.o : $(LIBDIR)/empirical.h \
	$(LIBDIR)/empirical.c \
	$(LIBDIR)/modeling.h
	$(CC) -c $(LIBDIR)/empirical.c -o $(BUILDDIR)/empirical.o

$(BUILDDIR)/nuclear_matter.o : $(LIBDIR)/nuclear_matter.h \
	$(LIBDIR)/nuclear_matter.c \
	$(LIBDIR)/mathconst.h \
	$(LIBDIR)/phyconst.h
	$(CC) -c $(LIBDIR)/nuclear_matter.c -o $(BUILDDIR)/nuclear_matter.o

$(BUILDDIR)/nuclear_surface_en.o : $(LIBDIR)/nuclear_surface_en.h \
	$(LIBDIR)/nuclear_surface_en.c \
	$(LIBDIR)/modeling.h
	$(CC) -c $(LIBDIR)/nuclear_surface_en.c -o $(BUILDDIR)/nuclear_surface_en.o

$(BUILDDIR)/coulomb.o : $(LIBDIR)/coulomb.h \
	$(LIBDIR)/coulomb.c \
	$(LIBDIR)/empirical.h \
	$(LIBDIR)/mathconst.h \
	$(LIBDIR)/phyconst.h
	$(CC) -c $(LIBDIR)/coulomb.c -o $(BUILDDIR)/coulomb.o

$(BUILDDIR)/nuclear_en.o : $(LIBDIR)/nuclear_en.h \
	$(LIBDIR)/nuclear_en.c
	$(CC) -c $(LIBDIR)/nuclear_en.c -o $(BUILDDIR)/nuclear_en.o

$(BUILDDIR)/crust.o : $(LIBDIR)/crust.h \
	$(LIBDIR)/crust.c \
	$(LIBDIR)/modeling.h
	$(CC) -c $(LIBDIR)/crust.c -o $(BUILDDIR)/crust.o

$(BUILDDIR)/core.o : $(LIBDIR)/core.h \
	$(LIBDIR)/core.c \
	$(LIBDIR)/modeling.h
	$(CC) -c $(LIBDIR)/core.c -o $(BUILDDIR)/core.o

clean:
	$(RM) $(OBJS) $(APP)
