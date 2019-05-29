#############################################
#      THM Make Include File          #
#############################################

# OpenBlas and Lapack libraries
BLASDIR := /usr/local/opt/openblas
BLASINC := $(BLASDIR)/include
BLASLIB := $(BLASDIR)/lib
LAPACKDIR := /usr/local/opt/lapack
LAPACKINC := $(LAPACKDIR)/include
LAPACKLIB := $(LAPACKDIR)/lib
# fortran compiler and flags
FC := gfortran
FFLAGS := -O2
FINC := -I $(BLASDIR) -I $(LAPACKDIR)
FLINK := -L $(BLASLIB) -L $(LAPACKLIB)