# -*-makefile-*-

#INCLUDES += 

# Do not disable nan! All other optimizations are fine, since exp/pow are used in expressions with low accuracy coefficients
MATHFLAGS = -funsafe-math-optimizations -fno-math-errno -fno-rounding-math -fno-trapping-math -fno-signed-zeros -freciprocal-math

ifneq (,$(findstring debug,$(MAKECMDGOALS)))
  FFLAGS = -MD -fno-omit-frame-pointer -L../.. -J../../obj -O0 -ggdb -cpp
else
  FFLAGS = -MD -fno-omit-frame-pointer -L../.. -J../../obj -O2 -ggdb -Ofast -cpp
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
FFLAGS += $(MATHFLAGS) -fdiagnostics-color=auto -fbounds-check

.SUFFIXES: $(SUFFIXES) .cpp .f90

vpath %.cpp src
vpath %.f90 src

CSOURCES = $(wildcard src/*.cpp)
FSOURCES = $(wildcard src/*.f90)
COBJFILES = $(patsubst %.cpp, obj/%.o, $(notdir $(CSOURCES)))
FOBJFILES = $(patsubst %.f90, obj/%.o, $(notdir $(FSOURCES)))
OBJFILES = $(FOBJFILES) $(COBJFILES)
#LDFLAGS += -L../../

LIBS +=	 -lcpr \
	-lcurl \
	-lstx-exparser \
	-lgfortran -lm -lgomp \
        -lpthread \
        -ljsoncpp \
        -lroadsurf	

all: roadrunner
debug: roadrunner
release: roadrunner

format:
	clang-format -i -style=file src/*.h src/*.cpp

install:
	mkdir -p $(bindir)
	$(INSTALL_PROG) roadrunner $(bindir)/

clean:
	rm -f obj/*.o obj/*.d *.mod roadrunner

roadrunner: objdir $(OBJFILES)
	$(CXX) -L../.. -o $@ $(OBJFILES) $(LIBS)

objdir:
	@mkdir -p obj

obj/%.o : %.cpp
	$(CXX) -c $(CFLAGS) -MD -MF $(patsubst obj/%.o, obj/%.d, $@) -MT $@ -o $@ $<

obj/SmartMetSource.o: FLAGS += -Wno-deprecated-declarations
obj/JsonTools.o: FLAGS += -Wno-deprecated-declarations

obj/%.o : %.f90
	ccache gfortran $(FFLAGS) -c -o $@ $< $(LIBS)

-include Dependencies

ifneq ($(wildcard obj/*.d),)
-include $(wildcard obj/*.d)
endif
