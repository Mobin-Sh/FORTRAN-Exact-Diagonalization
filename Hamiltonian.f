	Module Hamiltonian

	contains
		!----------------------------------------------------------------------------------------------------	  
		subroutine Cd(A,site)
			integer, parameter :: LKind = selected_int_kind (10) !Digits used in SN
			integer (kind=LKind) :: A,site
			
			A = A+site
			
		End subroutine Cd 
		!----------------------------------------------------------------------------------------------------	
		  
		subroutine C(A,site)
			integer, parameter :: LKind = selected_int_kind (10) !Digits used in SN
			integer (kind=LKind) :: A,site
			
			A = A-site
			
		End subroutine C
		!----------------------------------------------------------------------------------------------------		

		!----------------------------------------------------------------------------------------------------	
		subroutine Binarize (num,bin,n)
			integer, parameter :: LKind = selected_int_kind (10) !Digits used in SN
			integer (kind=LKind) :: num,dec !Left State Number, Right State Number, TwoExpontents
			Logical bin(n)
			Integer :: n
		
			dec=num

			DO i=1,n
				if (mod(dec,2)==0) then
					bin(i)=.False.
				else
					bin(i)=.True.
					
				end if
				dec=dec/2 !Notice the use of truncated result
				counter=counter+1
			
				if (dec==0) then
					exit
				end if
			End Do
			

		End subroutine Binarize 
		!----------------------------------------------------------------------------------------------------	        	

		!----------------------------------------------------------------------------------------------------------	
		
		
		subroutine HG(H,n)

		     !Implicit none
		     
		     !----------------------------------------------------------------------------------------
			!Variables  
			Integer :: n
			real :: t=3,g=4
			Integer i,j
			integer, parameter :: LKind = selected_int_kind (10) !Digits used in SN
			integer (kind=LKind) :: LSN,RSN,te(n) !Left State Number, Right State Number, TwoExpontents
			complex :: H(2**(n/2),2**(n/2))
			character(32) :: gchar
			Logical binSN(n)
			
			
			
			!----------------------------------------------------------------------------------------
			!Pre-Calculations
			
			!Defining Numbers 2**n Array
			Do i=1, n
				te(i)=2**(i-1)
			End do
			
			
			!Zeroing the Hamiltonian (Non-Sparse only)
			Do i=1, 2**(n/2)
				Do j=1, 2**(n/2)
					H(i,j)=0
				End do
			End do
			
			!Zeroing the Binary Format
			Do i=1, n
				binSN(i)=.False.
			End do
			

			!----------------------------------------------------------------------------------------
			!Defining The Hamiltonian
			
			Do RSN=0,(2**(n/2))-1
				!Binarize The RSN in binSN
				Call Binarize(RSN,binSN,n)
				Do j=1, n
				
			!Hamiltonian Terms:
			
			!on-Site
					LSN=RSN
					If(binSN(j) .eqv. .True.) then !checking wheter j is filled
						!call Cd(LSN,te(j))  !No need for that
						!call C(LSN,te(j))   !No need for that
						H(LSN+1,RSN+1)= H(LSN+1,RSN+1) + 2*t + 4*g
						
						!This is because of the Open Boudnary Condition of two term interaction of g
						!If (j==1 .or. j==n) then
						!	H(LSN+1,RSN+1)= H(LSN+1,RSN+1) - 2*g !OBC
						!End If
						
					End If
					
						
			!Skipping the Last sites == For the j,j+1 Terms
					If(j<n) then	
					
					!Hoppings 1
						LSN=RSN
						If( binSN(j+1)  .and. (.not. binSN(j))) then !checking wheter j is empty and j+1 is filled
							call Cd(LSN,te(j)) 
							call C(LSN,te(j+1)) 
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1)-t
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1)-t !Hermitian Conjugate
						End If	


					!Pairing Term 1
						LSN=RSN
						If( binSN(j+1)  .and. ( binSN(j))) then !checking wheter j and j+1 are filled
							call C(LSN,te(j)) 
							call C(LSN,te(j+1)) 
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) +t
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1) +t !Hermitian Conjugate
						End If
						
						
					!Interaction Term 1
						LSN=RSN
						If( binSN(j) .and. binSN(j+1)) then !checking wheter j and j+1 are filled
							!call C(LSN,te(j+1))   !No need for that
							!call Cd(LSN,te(j+1))  !No need for that
							!call C(LSN,te(j))     !No need for that
							!call Cd(LSN,te(j))    !No need for that
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) - 4*g
						End If		
						
					End If				


			!Skipping the Last 2 sites== For the j,j+2 Terms		
					If(j<n-1) then			
					
					
					!Hoppings 2
						LSN=RSN
						If( binSN(j+2)  .and. (.not. binSN(j))) then !checking wheter j is empty and j+2 is filled
							call Cd(LSN,te(j)) 
							call C(LSN,te(j+2)) 
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1)-g
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1)-g !Hermitian Conjugate
						End If
						
						
					!Pairing Term 2-a
						LSN=RSN
						If( binSN(j+2)  .and. ( binSN(j))) then !checking wheter j and j+2 are filled
							call C(LSN,te(j)) 
							call C(LSN,te(j+2)) 
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) +g
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1) +g !Hermitian Conjugate
						End If	
						
						
					!Pairing Term 2-b
						LSN=RSN
						If( binSN(j+1)  .and. ( binSN(j))) then !checking wheter j and j+1 are filled
							call C(LSN,te(j)) 
							call C(LSN,te(j+1)) 
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) +t
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1) +t !Hermitian Conjugate
						End If	
						
						
					!Interaction Term 2-a
						LSN=RSN
						If( binSN(j+2) .and. binSN(j+1) .and. binSN(j)) then !checking wheter j and j+1 are filled
							call C(LSN,te(j+2))  
							!call C(LSN,te(j+1))  !No need for that
							!call Cd(LSN,te(j+1))     !No need for that
							call C(LSN,te(j))
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) - g
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1) - g !Hermitian Conjugate
						End If		
						
						
					!Interaction Term 2-b
						LSN=RSN
						If( (.not. binSN(j)) .and. binSN(j+1) .and. binSN(j+2)) then !checking wheter j and j+1 are filled
							call C(LSN,te(j+2))
							!call C(LSN,te(j+1))  !No need for that
							!call Cd(LSN,te(j+1))     !No need for that
							call Cd(LSN,te(j))
							H(LSN+1,RSN+1)= H(LSN+1,RSN+1) + g
							H(RSN+1,LSN+1)= H(RSN+1,LSN+1) + g !Hermitian Conjugate
						End If		
						
					End If				
					
					
					!PBC
					
					

					
					!	LSN=RSN 
					!	If( binSN(j) .and. (.not. binSN(j+1))) then
					!		call Cd(LSN,te(j+1)) 
					!		call C(LSN,te(j)) 
					!		H(LSN+1,RSN+1)= H(LSN+1,RSN+1)-t
					!	End If
					
					
				End do
			End do	
			
			Print *,"Hamiltonian Made Successfully!"
			!Call MatrixForm(H,2**n)
			
			
			
	  

		End subroutine
		
	End Module Hamiltonian	
