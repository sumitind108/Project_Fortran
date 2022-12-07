project_1: H_1D_serial.f90 in_1.dat
	gfortran H_1D_serial.f90 -o project_1
	./project_1 < in_1.dat > out_1.dat

project_2: H_2D_serial.f90 in_2.dat
	gfortran H_2D_serial.f90 -o project_2
	./project_2 < in_2.dat > out_2.dat

project_3: H_3D_serial.f90 in_3.dat
	gfortran H_3D_serial.f90 -o project_3
	./project_3 < in_3.dat > out_3.dat

clean:
	rm -rf out*.dat project*