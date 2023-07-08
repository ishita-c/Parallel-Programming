compile:
	nvcc -Wno-deprecated-gpu-targets -std=c++11 -O3 -gencode=arch=compute_35,code=sm_35 main.cu -o exec
run:
	./exec test_case/input1.bin test_case/input2.bin output.bin
large:
	./exec test_case_large/input1.bin test_case_large/input2.bin output_large_9.bin
check:
	python3 checker.py -f out_gold_large.bin output_large_9.bin -e 4
clean:
	rm exec
