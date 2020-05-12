	Module EVPkg
	
	!Implicit none

	contains
!----------------------------------------------------------------------------------------------------	     
		subroutine MatrixForm(MFMatrix,MFDimension)
		!Subroutine to print the matrix in a "nice" matrix form.
		!Entry is : (MFMatrix=Matrix you want to print,MFDimension=its' Dimension)
			complex :: MFMatrix(MFDimension,MFDimension)
			Write(*,*)"Matrix Form:"
			Do, i=1,MFDimension
				Do, j=1, MFDimension
					write(*, '("	",E22.16,",",E22.16)', advance="no") MFMatrix(i,j)
				End do
				write(*, *) ""
			End do
		
		End subroutine MatrixForm 
!----------------------------------------------------------------------------------------------------	
		
		subroutine AllEV(A, M, EVL)
		!A is the Matrix | M is the order | EVList is an array of all of the EVs, Output
	
			!Variable Declaration
			Integer :: M
			Complex AA(M,M), A(M,M), EVL(M), b(M), DUMMY(1,1), WORK(2*M)
		
			AA=A
			
			Do i=1, M
				Print *, AA(i,M)
			End do
			
			!Using Lapack EV calculator
			call ZGEEV('N', 'N', M, AA, M, b, DUMMY, 1, DUMMY, 1, WORK, 2*M, WORK, ok)
		
			Do i=1, M
				Print *, "EV",i,"	",b(i)
			End do
			

			
			!Printing The Result
			If (ok .eq. 0) then
				Do i=1, M
					EVL(i)= b(i)
				End do
				Write (*,*) "EVList was successfully calculated."
			Else
				Write (*,*) "An error occured, Please check your Entry Matrix."
			Endif
			Print *, "Am Here."
		END subroutine AllEV
!----------------------------------------------------------------------------------------------------		
		
	End Module EVPkg
