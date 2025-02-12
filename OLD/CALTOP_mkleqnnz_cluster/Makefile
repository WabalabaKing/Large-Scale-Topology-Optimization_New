# Define paths for SPOOLES and ARPACK
SPOOLES_PATH = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/OLD/Deps/Spooles/src
ARPACK_PATH = /home/prateek/Documents/Github/Large-Scale-Topology-Optimization/OLD/Deps/ARPACK

# Compiler and flags
CFLAGS = -Wall -O3 -I $(SPOOLES_PATH) -DARCH="Linux" -DSPOOLES -DARPACK -DMATRIXSTORAGE -DNETWORKOUT
FFLAGS = -Wall -O3 -fallow-argument-mismatch

CC = cc
FC = gfortran

#Directories
SRC_DIR = src


# Object file rules
.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

# Include additional makefile if necessary
include Makefile.inc

# Source files and objects
SCCXMAIN = $(SRC_DIR)/ccx_2.15.c

OCCXF = $(SCCXF:$(SRC_DIR)/%.f=%.o)
OCCXC = $(SCCXC:$(SRC_DIR)/%.c=%.o)
#OCCXMAIN = $(SCCXMAIN:.c=.o)
OCCXMAIN = ccx_2.15.o

#OCCXF = $(SCCXF:.f=.o)
#OCCXC = $(SCCXC:.c=.o)
#OCCXMAIN = $(SCCXMAIN:.c=.o)

# Libraries
LIBS = \
	$(SPOOLES_PATH)/spooles.a \
	$(ARPACK_PATH)/libarpack_INTEL.a \
	-lpthread -lm -lc


# Object file compilation rules
%.o: $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -c $< -o $(notdir $@)

# Build executable
calTop.exe: $(OCCXMAIN) ccx_2.15.a $(LIBS)
	./date.pl; $(CC) $(CFLAGS) -c $(SRC_DIR)/ccx_2.15.c; $(FC) -Wall -O3 -o $@ $(OCCXMAIN) ccx_2.15.a $(LIBS)
	@echo "Build successful!"

# Build static library
ccx_2.15.a: $(OCCXF) $(OCCXC)
	ar vr $@ $?


# Clean up binaries and object files
clean:
	rm -f *.o calTop.exe ccx_2.15.a
