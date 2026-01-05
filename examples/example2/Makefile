# -*-makefile-*-

REQUIRES = gdal jsoncpp libpqxx fmt

include $(shell echo $${PREFIX-/usr})/share/smartmet/devel/makefile.inc

INCLUDES += -I/usr/include/smartmet

# Do not disable nan! All other optimizations are fine, since exp/pow are used in expressions with low accuracy coefficients
MATHFLAGS = -funsafe-math-optimizations -fno-math-errno -fno-rounding-math -fno-trapping-math -fno-signed-zeros -freciprocal-math

FFLAGS = -I/usr/include/roadsurf -Isrc -Jobj -MD -fno-omit-frame-pointer -ggdb -cpp

ifneq (,$(findstring debug,$(MAKECMDGOALS)))
  FFLAGS += -O0
else
  FFLAGS += -O2 -Ofast
endif

ifeq ($(ASAN), yes)
  FFLAGS += -fsanitize=address -fsanitize=pointer-compare -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize-address-use-after-scope
endif

ifeq ($(TSAN), yes)
  FFLAGS += -fsanitize=thread
endif

ifeq ($(BSAN), yes)
  FFLAGS += -fbounds-check
endif

CFLAGS += $(MATHFLAGS) #-D_GLIBCXX_DEBUG
FFLAGS += $(MATHFLAGS) -fdiagnostics-color=$(GCC_DIAG_COLOR) -fbounds-check

.SUFFIXES: $(SUFFIXES) .cpp .f90

vpath %.cpp src
vpath %.f90 src

CSOURCES = $(wildcard src/*.cpp)
FSOURCES = $(wildcard src/*.f90)
COBJFILES = $(patsubst %.cpp, obj/%.o, $(notdir $(CSOURCES)))
FOBJFILES = $(patsubst %.f90, obj/%.o, $(notdir $(FSOURCES)))
OBJFILES = $(FOBJFILES) $(COBJFILES)

LIBS +=	-L$(libdir) \
	-lsmartmet-macgyver \
	-lsmartmet-newbase \
	-lsmartmet-locus \
	-lroadsurf \
	-lcpr \
	-lcurl \
	-lstx-exparser \
	-lboost_program_options \
	-lboost_filesystem \
	-lboost_date_time \
	-lboost_iostreams \
	-lboost_thread \
	-lboost_timer \
	$(REQUIRED_LIBS) \
	-lgfortran -lm -lgomp \
	-lpthread \
        -lboost_serialization

all: roadrunner
debug: roadrunner
release: roadrunner

rpm: clean smartmet-roadrunner.spec
	rm -f smartmet-roadrunner.tar.gz # Clean a possible leftover from previous attempt
	tar -czvf smartmet-roadrunner.tar.gz --exclude test --exclude-vcs --transform "s,^,smartmet-roadrunner/," *
	rpmbuild -tb smartmet-roadrunner.tar.gz
	rm -f smartmet-roadrunner.tar.gz

format:
	clang-format -i -style=file src/*.h src/*.cpp

install:
	mkdir -p $(bindir)
	$(INSTALL_PROG) roadrunner $(bindir)/

clean:
	rm -f obj/*.o obj/*.d *.mod roadrunner

roadrunner: objdir $(OBJFILES)
	$(CXX) $(LDFLAGS) -o $@ $(OBJFILES) $(LIBS)

objdir:
	@mkdir -p obj

obj/%.o : %.cpp
	$(CXX) -c $(CFLAGS) $(INCLUDES) -MD -MF $(patsubst obj/%.o, obj/%.d, $@) -MT $@ -o $@ $<

obj/SmartMetSource.o: FLAGS += -Wno-deprecated-declarations
obj/JsonTools.o: FLAGS += -Wno-deprecated-declarations

obj/%.o : %.f90
	ccache gfortran $(FFLAGS) -c -o $@ $<

-include Dependencies

ifneq ($(wildcard obj/*.d),)
-include $(wildcard obj/*.d)
endif
