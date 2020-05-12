all: EVPkg.f
	mkdir ./Output/ -p
	mkdir ./misc/ -p
	gfortran -c EVPkg.f -llapack -ffixed-line-length-100
	gfortran -c Hamiltonian.f 
	gfortran Run.f EVPkg.o Hamiltonian.o -o./Output/output.o -llapack -ffixed-line-length-100 
	cp ./evpkg.mod	./misc/
	cp ./EVPkg.o	./misc/
	rm evpkg.mod EVPkg.o
	./Output/output.o

