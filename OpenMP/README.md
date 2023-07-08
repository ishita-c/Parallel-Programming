Number of approaches: 3

Approach 1
1:Idea: Brute Force, O(n^3) based algorithm with parallel read and write
1:Why attempted: I attempted this approach to get a baseline for further optimizations. Reading from input binary 
file and writing to output binary file was done parallely using tasks equal to number of threads. Contiguous chunks 
of binary file were allocated to each thread thereby preventing false sharing. The algorithm for matrix multiplication
used was O(n^3) without actually storing the entire matrix in memory. 
1:Results (Specify speedup): The results with this approach were very poor and the speedup was negligible because O(n^3)
algorithm is very computationally expensive and not very practical thereby nullifying the parallelization benefits from 
reading and writing. 
1:Drawback: The major drawback of this approach was the O(n^3) algorithm for calculation of the resultant matrix.
Also, the algorithm doesn't use the fact that the matrices are sparse and thus loses out on memory and runtime 
optimization opportunities.


Approach 2
2:Idea: Parallelized O(k^2) based algorithm with parallel read and write. In this algorithm, for every non-0 input 
block, we iterate over all non-0 input blocks and keep aggregating the contributions. Contiguous sets of rows are 
assigned to each thread for computation to avoid false sharing. This also allows for easy splitting of the total 
work into disjoint tasks.
2:Why attempted: I made the observation that only pairs of non-0 input blocks can make some contribution to the final 
resultant matrix and thus attempted a O(k^2) based algorithm for computation divided into tasks equal to the number 
of threads. 
2:Results (Specify speedup): The results were fairly good using this approach. The overall runtime for the large test 
case provided with about 20k non-0 input blocks came down to 50s with 1 core and 18s with 16 cores, thereby giving a 
speedup of upto 3X.
2:Drawback: The major drawback here is that for every non-0 input block, only some of the other blocks are compatible 
for contributing to the final output and thus it is wasteful to iterate over all the non-0 blocks.

Final Approach
3:Idea: Parallelized O(k*d) based algorithm with parallel read and write. In this approach we do some pre-computation 
such as storing the non-0 input blocks by their row index in an unordered map. Thus for every non-0 input block, 
there will be atmost d = n/m non-0 blocks to iterate over. Operations while are polynomial time in "m" are considered 
to be constant in this complexity analysis. Reading, calculation and writing, all the three use tasks equal to number 
of threads, with work divided into contiguous blocks of memory thus ensuring that the program is parallel and scalable.
3:Why attempted: This approach is the most complexity efficient among all the approaches explored so far and thus is 
the most suitable candidate for my final solution.
3:Results (Specify speedup): The runtime results are much better than the previous approaches.The overall runtime for
the large test case provided with about 20k non-0 input blocks came down to 20s with 1 core and 3.75s with 16 cores,
thereby giving a speedup of upto 5X. For even larger inputs such as 65k Non-0 input blocks, the speedup was upto 8X
with 16 cores. Thus we observe that this final approach is fairly well parallel and scalable.
3:Drawback: One minor drawback of this approach is that the parallel execution may not be utilized to its full 
capacity if the data in the sparse distribution is not random and concentrated in one part of the matrix. 

--------------------------------------------------------------------------------------------------------------
Final scalability analysis (time in milliseconds)

Non-0 input blocks  |   Non-0 output blocks  |   1 core    |   2 cores  |   4 cores    |    8 cores   |   16 cores    |
    2^10 =  1024    |          2119          |      242    |      139   |       83     |       75     |       99      |       
    2^11 =  2048    |          4749          |      579    |      353   |      199     |      158     |      201      |              
    2^12 =  4096    |         11949          |     1494    |      974   |      555     |      346     |      423      |
    2^13 =  8192    |         37230          |     4418    |     3159   |     1783     |      973     |      981      |
    2^14 = 16384    |        136647          |    14851    |    10933   |     6424     |     3537     |     2455      |
    2^15 = 32768    |        521526          |    54033    |    40309   |    23832     |    13528     |     7035      |
    2^16 = 65536    |       1884747          |   202716    |   155697   |    92557     |    49646     |    26241      |


Comments: We observe that for smaller inputs(less that 2^13 Non-0 input blocks), parallel execution with 8 threads 
gives the lowest runtime whereas for larger inputs, execution with 16 threads is faster. This is probably because 
the overhead of synchronization more than compensates for the benefits of parallelization while switching from 8 to 
16 threads for smallers input sizes. For larger inputs,this overhead is not significant against the effects of 
parallelization. We observe upto 8X speedup in case of large test case with 2^16 Non-0 input blocks and expect this 
speedup to increase further for even larger inputs.

--------------------------------------------------------------------------------------------------------------
