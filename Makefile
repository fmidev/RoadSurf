# -*-makefile-*-

# Installation

INSTALL_PROG = install -p -m 775
INSTALL_DATA = install -p -m 664

LDFLAGS += -shared -rdynamic

ifeq ($(origin PREFIX), undefined)
  PREFIX = /usr
else
  PREFIX = $(PREFIX)
endif

processor := $(shell uname -p)
ifeq ($(processor), x86_64)
  libdir ?= $(PREFIX)/lib64
else
  libdir ?= $(PREFIX)/lib
endif

includedir ?= $(PREFIX)/include

objdir = obj

# Do not disable nan! All other optimizations are fine, since exp/pow are used in expressions with low accuracy coefficients
MATHFLAGS = -funsafe-math-optimizations -fno-math-errno -fno-rounding-math -fno-trapping-math -fno-signed-zeros -freciprocal-math

FC = gfortran

FLAGS = -MD -fno-omit-frame-pointer -Jobj -ggdb -cpp -fPIC

ifneq (,$(findstring debug,$(MAKECMDGOALS)))
  FLAGS += -O0 -Wall
else
  FLAGS += -O2 -Ofast
endif

ifeq ($(ASAN), yes)
  FLAGS += -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize-address-use-after-scope
endif

ifeq ($(TSAN), yes)
  FLAGS += -fsanitize=thread
endif

ifeq ($(BSAN), yes)
  FLAGS += -fbounds-check
endif

FLAGS += $(MATHFLAGS) -fdiagnostics-color=auto

.SUFFIXES: $(SUFFIXES) .f90
.PHONY: all debug release clean

vpath %.f90 src
vpath %.o obj
vpath %.mod obj

LIBFILE = libroadsurf.so
LIBS +=	-L$(libdir) -l$(FC)  
SRCS = $(wildcard src/*.f90)
OBJS = $(patsubst %.f90, %.o, $(notdir $(SRCS)))

all: objdir $(LIBFILE)
debug: all
release: all

$(LIBFILE): $(OBJS:%.o=$(objdir)/%.o) | objdir
	$(FC) $(LDFLAGS) $^ $(LIBS) -o $(LIBFILE)

$(objdir)/%.o: %.f90
	$(FC) $(FLAGS) -c $< -o $@

$(objdir)/RoadSurf.o: $(objdir)/RoadSurfVariables.o
$(objdir)/BalanceModel.o: $(objdir)/RoadSurf.o
$(objdir)/BoundaryLayer.o: $(objdir)/RoadSurf.o
$(objdir)/Cond.o: $(objdir)/RoadSurf.o
$(objdir)/ConnectFortran2Carrays.o: $(objdir)/RoadSurf.o
$(objdir)/Coupling.o: $(objdir)/RoadSurf.o
$(objdir)/Initialization.o: $(objdir)/RoadSurf.o
$(objdir)/InputOutput.o: $(objdir)/RoadSurf.o
$(objdir)/ModRadiation.o: $(objdir)/RoadSurf.o
$(objdir)/Relaxation.o: $(objdir)/RoadSurf.o
$(objdir)/Storage.o: $(objdir)/RoadSurf.o
$(objdir)/SunPosition.o: $(objdir)/RoadSurf.o

objdir:
	@mkdir -p $(objdir)

examples: example1 example2

example1: all
	$(MAKE) -C examples/example1

example2: all
	$(MAKE) -C examples/example2

rpm : clean roadsurf.spec
	rm -f roadsurf.tar.gz # Clean a possible leftover from previous attempt
	tar -czvf roadsurf.tar.gz --exclude test --exclude examples --exclude-vcs --transform "s,^,roadsurf/," *
	rpmbuild -tb roadsurf.tar.gz
	rm -f roadsurf.tar.gz

clean:
	rm -rf $(LIBFILE) $(objdir)
	$(MAKE) -C examples/example1 clean
	$(MAKE) -C examples/example2 clean

install:
	mkdir -p $(libdir)
	$(INSTALL_PROG) $(LIBFILE) $(libdir)/$(LIBFILE)
	mkdir -p $(includedir)/roadsurf
	$(INSTALL_DATA) $(objdir)/*.mod $(includedir)/roadsurf/
	$(INSTALL_DATA) src/Constants.h $(includedir)/roadsurf/
