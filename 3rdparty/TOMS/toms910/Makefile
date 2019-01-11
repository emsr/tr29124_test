###############################################################
#
# Makefile for e_float
#
###############################################################


# ------------------------------------------------------------------------------
# filelist
# ------------------------------------------------------------------------------
include MakefileFiles.gmk

ifeq ($(MP),)
  TYPE_MP = efx
endif

ifeq ($(MP),efx)
  TYPE_MP = efx
endif

ifeq ($(MP),f90)
  TYPE_MP = f90
endif

ifeq ($(MP),gmp)
  TYPE_MP = gmp
endif

ifeq ($(MP),mpfr)
  TYPE_MP = mpfr
endif

ifeq ($(MP),pyd)
  TYPE_MP = pyd
endif

# ------------------------------------------------------------------------------
# symbols
# ------------------------------------------------------------------------------
PAR_OPEN  = \(
PAR_CLOSE = \)


# ------------------------------------------------------------------------------
# targets
# ------------------------------------------------------------------------------
PATH_TARGET  = unix-$(TYPE_MP)
PATH_OBJ     = $(PATH_TARGET)/obj
PATH_ERR     = $(PATH_TARGET)/err
TARGET       = e_float
EXTENSION    = exe


# ------------------------------------------------------------------------------
# Tool parameters
# ------------------------------------------------------------------------------
#CFLAGS_BASE = -std=c++0x -Wall -ffast-math -I$(CURDIR)
CFLAGS_BASE = -Wall -ffast-math -I$(CURDIR)
LFLAGS_BASE =

ifeq ($(TYPE_MP),efx)
  CFLAGS_BASE := $(CFLAGS_BASE) -DE_FLOAT_TYPE_EFX
  FILES_MP = e_float/efx/e_float_efx
endif

ifeq ($(TYPE_MP),f90)
  CFLAGS_BASE := $(CFLAGS_BASE) -DE_FLOAT_TYPE_F90
  FILES_MP = e_float/f90/e_float_f90
endif

ifeq ($(TYPE_MP),gmp)
  CFLAGS_BASE := $(CFLAGS_BASE) -DE_FLOAT_TYPE_GMP
  LFLAGS_BASE := $(LFLAGS_BASE) -lgmp
  FILES_MP = e_float/gmp/e_float_gmp
endif

ifeq ($(TYPE_MP),mpfr)
  CFLAGS_BASE := $(CFLAGS_BASE) -DE_FLOAT_TYPE_MPFR
  LFLAGS_BASE := $(LFLAGS_BASE) -lgmp -lmpfr
  FILES_MP = e_float/mpfr/e_float_mpfr
endif

ifeq ($(TYPE_MP),pyd)
  CFLAGS_BASE := $(CFLAGS_BASE) -fPIC -I/usr/include/boost -I/usr/include/python2.6 -DE_FLOAT_TYPE_MPFR
  LFLAGS_BASE := $(LFLAGS_BASE) -shared -lgmp -lmpfr -lboost_python
  TARGET      := e_float_pyd
  EXTENSION   := so
  FILES_MP         = e_float/mpfr/e_float_mpfr
  FILES_INTEROP   := interop/pyd/e_float_pyd
  FILES_TEST      :=
  FILES_TEST_SPOT :=
  FILES_TEST_REAL :=
  FILES_TEST_IMAG :=
  FILES_EXAMPLES  :=
endif


CFLAGS = -O3 $(CFLAGS_BASE)
LFLAGS = $(LFLAGS_BASE)


# ------------------------------------------------------------------------------
# object files
# ------------------------------------------------------------------------------
FILES_ALL = $(FILES_MP)        \
            $(FILES_EF)        \
            $(FILES_UTILITY)   \
            $(FILES_FUNCTIONS) \
            $(FILES_EXAMPLES)  \
            $(FILES_INTEROP)   \
            $(FILES_TEST)      \
            $(FILES_TEST_SPOT) \
            $(FILES_TEST_REAL) \
            $(FILES_TEST_IMAG)

FILES_O = $(addprefix $(PATH_OBJ)/, $(notdir $(addsuffix .o, $(FILES_ALL))))


# ------------------------------------------------------------------------------
# VPATH definition: VPATH is required for make to find the source files.
# ------------------------------------------------------------------------------
VPATH    := $(sort $(dir $(FILES_ALL)))


# ------------------------------------------------------------------------------
# Development tools
# ------------------------------------------------------------------------------
GCC   = gcc
CC    = g++
SED   = sed
ECHO  = echo
MAKE  = make
MKDIR = mkdir
RM    = rm
GAWK  = gawk


# ------------------------------------------------------------------------------
# Main-Dependencies (target: all)
# ------------------------------------------------------------------------------
.PHONY: all
all: OutputDirs Version $(PATH_TARGET)/$(TARGET).$(EXTENSION)


# ------------------------------------------------------------------------------
# Main-Dependencies (target: rebuild)
# ------------------------------------------------------------------------------
.PHONY: rebuild
rebuild: clean all


# ------------------------------------------------------------------------------
# Clean project (target: clean)
# ------------------------------------------------------------------------------
.PHONY: clean
clean: Version
	@$(ECHO) +++ cleaning all
	@-$(RM) -f $(PATH_TARGET)/$(TARGET).$(EXTENSION) > NUL  2> NUL
	@-$(RM) -f $(PATH_OBJ)/*.*                       > NUL  2> NUL
	@-$(RM) -f $(PATH_ERR)/*.*                       > NUL  2> NUL


# ------------------------------------------------------------------------------
# create output directories
# ------------------------------------------------------------------------------
.PHONY: OutputDirs
OutputDirs:
  # Create the output directories
	@-$(MKDIR) -p $(PATH_TARGET)
	@-$(MKDIR) -p $(PATH_OBJ)
	@-$(MKDIR) -p $(PATH_ERR)


# ------------------------------------------------------------------------------
# version information
# ------------------------------------------------------------------------------
.PHONY: Version
Version:
  # Print the GNU make version and the GCC version
	@-$(MAKE) --version | $(GAWK) '{ if(match($$0, "Make") > 0) print("GNUmake version: " $$0) }'
	@-$(GCC)  --version | $(GAWK) '{ if(match($$0, "gcc")  > 0) print("GCC version: " $$0) }'


# ------------------------------------------------------------------------------
# link/locate application
# ------------------------------------------------------------------------------
$(PATH_TARGET)/$(TARGET).$(EXTENSION) : $(FILES_O)
	@$(ECHO) +++ linking application: link to $(PATH_TARGET)/$(TARGET).$(EXTENSION)
	@$(CC) $(FILES_O) $(LFLAGS) -o $(PATH_TARGET)/$(TARGET).$(EXTENSION)


# ------------------------------------------------------------------------------
# include dependency files
# ------------------------------------------------------------------------------
ifneq ($(MAKECMDGOALS),rebuild)
-include $(subst .o,.d, $(FILES_O))
endif


# ------------------------------------------------------------------------------
# pattern rule for cpp-files
# ------------------------------------------------------------------------------
$(PATH_OBJ)/%.o : %.cpp
  # Compile the source file...
  # ...and Reformat (using sed) any possible error/warning messages for the Visual Studio build window
  # ...and Generate a dependency file (using the g++ flag -MM), reformating it with sed
	@-$(ECHO) +++ compiling: $<
	@-$(CC) $(CFLAGS) $< -c -o $(PATH_OBJ)/$(basename $(@F)).o 2> $(PATH_ERR)/$(basename $(@F)).err
	@-$(SED) -e 's|.h:\([0-9]*\),|.h(\1) :|' -e 's|:\([0-9]*\):|(\1) :|' $(PATH_ERR)/$(basename $(@F)).err
	@-$(CC) $(CFLAGS) -MM $< | $(SED) -e 's|\($(basename $(@F)).o\):|$(PATH_OBJ)/\1:|' > $(PATH_OBJ)/$(basename $(@F)).d
