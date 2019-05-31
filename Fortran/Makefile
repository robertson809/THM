#############################################
#           THM Make File Root directory   #
#############################################
include make.inc

install: tmp_eigsolver
	
tmp_eigsolver:
	$(FC) $(FFLAGS) -o tmp_driver fpml.f90 tmp_eigsolver.f90 tmp_driver.f90 $(FINC) $(FLINK) -lopenblas -llapack
	
tridiag:
	$(FC) $(FFLAGS) -o tridiag_driver tridiag.f90 tridiag_driver.f90
	
uninstall: clean
	@rm -f tridiag_driver
	@rm -f tmp_driver
	
clean:
	@rm -f *.mod