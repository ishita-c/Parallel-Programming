compile:
	mpic++ -std=c++17 -fopenmp -Ofast main.cpp -o a3
run:
	mpirun -np 8 ./a3 --taskid=2 --inputpath=test_cases/test3/test-input-3.gra --headerpath=test_cases/test3/test-header-3.dat --outputpath=test_cases/test3/out2_verbose.txt --verbose=1 --startk=3 --endk=3 --p=4
clean:
	rm a3
