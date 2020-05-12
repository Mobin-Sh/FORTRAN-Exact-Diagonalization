!for a single RUN: gfortran Run.f -llapack -ffixed-line-length-100 (c/p the subroutines from module)
!File EVPkg.f should be with this file as the source and module.
	  
	Program Main
	
          use EVPkg
          use Hamiltonian
		Implicit none
          
          !Variables  
		Integer,parameter :: n=6
		integer,parameter :: powern= 2**(n/2)
		Integer i
		complex :: H(powern,powern), EVlist(powern)
		
		
		!DATA
		Call HG(H,n)
		
		
		!Calling the Calculations
		Call MatrixForm(H,powern) !Print The Matrix Form of the Hamiltonian
		Call AllEV(H, powern, EVlist) !Calculate the Hamiltonian and puts it in EVlist as an array
		
		Print *, "Am Here."
		
		Do i=1, 2
			Print *, "EV",i,"	",EVlist(i)
		End do
  

	End program Main
