LD_LIBRARY_PATH=.
export LD_LIBRARY_PATH

all:
	g++ -std=c++11 -fopenmp -Wall -Ofast main.cpp -o exec -L. -lrary
run:
	./exec input my_out
clean:
	rm exec 