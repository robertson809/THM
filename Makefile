#############################################
#           THM Make File Root directory   #
#############################################
include make.inc

install: tridiag
	
tridiag:
	$(FC) $(FFLAGS) -o tridiag_driver tridiag.f90 tridiag_driver.f90
	
uninstall: clean
	@rm -f tridiag_driver
	
clean:
	@rm -f *.mod