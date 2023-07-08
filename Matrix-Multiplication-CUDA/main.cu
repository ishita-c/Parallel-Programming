#include <map>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace std::chrono;

int readFourBytes(fstream& fs){
    unsigned char a[4];
    fs.read((char *)&a, sizeof(a));
    int x = (int)a[0] | (int)a[1]<<8 | (int)a[2]<<16 | (int)a[3]<<24; 
    return x;
}

int readTwoBytes(fstream& fs){
    unsigned char a[2];
    fs.read((char *)&a, sizeof(a));
    int x = (int)a[0] | (int)a[1]<<8;
    return x;
}

void writeFourBytes(uint32_t x, fstream& fs){
    char* buffer;
    buffer = (char*)(&x);

    for(int i = 0; i < 4; i++)
        fs.write(buffer + i, sizeof(char));
}

void readData(int m, int d, int k, int *blocks, int* indices, int *mat, string input_file){
    int b = m*m;
    fstream fs(input_file, ios::in | ios::binary);
    fs.seekg(12, std::ios::beg);
    vector<pair<int, vector<int>>> data(k);
    for(int i_ = 0 ; i_ < k ; i_++){
        int i = readFourBytes(fs);
        int j = readFourBytes(fs);
        vector<int> temp(m*m);
        for(int k_ = 0; k_ < b; k_++){
            int x = readTwoBytes(fs);
            // mat[i_*b + k_] = x;
            temp[k_] = x;
        }
        // blocks[i_] = i*d + j;
        data[i_] = {i*d + j, temp};
    }
    sort(data.begin(), data.end());
    map<int, int> mp;
    for(int i=0;i<k;i++){
        blocks[i] = data[i].first;
        if(mp[blocks[i]/d] == 0){
            mp[blocks[i]/d] = 1;
            indices[blocks[i]/d] = i;
        }
        for(int j=0;j<b;j++){
            mat[i*b + j] = data[i].second[j];
        }
    }

    fs.close();
}

//TODO: Complete this
__global__
void calculateResult(int n, int m, int k1, int k2, int* blocks1, int* indices1, int* mat1, int* blocks2, int* indices2, int* mat2, int* blocksOutput, uint32_t* matOutput){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int d = n/m;
    int b = m*m;
    if(idx >= d*d) return;

    int x = idx/d;
    int y = idx%d;
    uint32_t MAX_VAL = pow(2,32) - 1;

    for(int i=indices1[x];i<k1;i++){
        int x1 = blocks1[i]/d;
        if(x1!=x) break;
        int y1 = blocks1[i]%d;
        for(int j=indices2[y1];j<k2;j++){
            int x2 = blocks2[j]/d;
            if(x2!=y1) break;
            int y2 = blocks2[j]%d;
            if(y2!=y) continue;

            int offset1 = i*b;
            int offset2 = j*b;

            for(int i_ = 0;i_<m;i_++){
                for(int j_ = 0;j_<m;j_++){
                    long val = 0;
                    for(int k=0;k<m;k++){
                        val += mat1[offset1 + i_*m + k]*mat2[offset2 + k*m + j_];
                        if (val >= (long)MAX_VAL){
                            val = (long)MAX_VAL;
                            break;
                        } 
                    }
                    if(val > 0){
                        long new_val = matOutput[idx*b + i_*m + j_] + val;
                        if (new_val >= (long)MAX_VAL){
                            new_val = (long)MAX_VAL;
                        }
                        matOutput[idx*b + i_*m + j_] = (uint32_t)new_val;
                        blocksOutput[idx] = 1;
                    }
                }
            }
        }
    }
}

// Outputter
void outputResult(int n, int m, int* blocksOutput, uint32_t* matOutput, string output_file){
    int d = n/m;
    int k = 0;
    int b = m*m;
    for(int idx=0;idx<d*d;idx++){
        if(blocksOutput[idx]) k++;
    }
    fstream fs(output_file, ios::out | ios::binary);
    writeFourBytes(n, fs);
    writeFourBytes(m, fs);
    writeFourBytes(k, fs);

    cout << "Output k = " << k << endl;

    for(int idx=0;idx<d*d;idx++){
        if(blocksOutput[idx] == 0) continue;
        writeFourBytes(idx/d, fs);
        writeFourBytes(idx%d, fs);
        // cout << idx/d << " " << idx%d << endl;
        for(int i=0;i<m*m;i++){
            writeFourBytes(matOutput[idx*b + i], fs);
            // cout << matOutput[idx*b + i] << " " ;
        }
        // cout << endl;
    }
    fs.close();
}

int main(int argc, char* argv[]){

    auto start = high_resolution_clock::now();

    string input_file1 = argv[1];
    string input_file2 = argv[2];
    string output_file = argv[3];
    
    fstream fs1(input_file1, ios::in | ios::binary);
    int n1 = readFourBytes(fs1);
    int m1 = readFourBytes(fs1);
    int k1 = readFourBytes(fs1);
    fs1.close();

    fstream fs2(input_file2, ios::in | ios::binary);
    int n2 = readFourBytes(fs2);
    int m2 = readFourBytes(fs2);
    int k2 = readFourBytes(fs2);
    fs2.close();

    cout << "k1" << " = " << k1 << endl;
    cout << "k2" << " = " << k2 << endl;

    if(n1 != n2 || m1 != m2){
        cout << "Error: Matrices have different dimensions!" << endl;
    }

    int n = n1;
    int m = m1;

    cout << "n" << " = " << n << endl;
    cout << "m" << " = " << m << endl;
    // dimension of block matrix
    int d = n/m;

    vector<pair<int, vector<int>>> data1(k1), data2(k2);

    int *blocks1;
    int *blocks2;
    int *indices1;
    int *indices2;
    int *mat1;
    int *mat2;
    int* blocksOutput;
    uint32_t* matOutput;

    cudaMallocManaged(&blocks1, k1*sizeof(int));
    cudaMallocManaged(&blocks2, k2*sizeof(int));
    cudaMallocManaged(&indices1, n*sizeof(int));
    cudaMallocManaged(&indices2, n*sizeof(int));
    cudaMallocManaged(&mat1, k1*m*m*sizeof(int));
    cudaMallocManaged(&mat2, k2*m*m*sizeof(int));

    cudaMallocManaged(&blocksOutput, d*d*sizeof(int));
    cudaMallocManaged(&matOutput, n*n*sizeof(uint32_t));

    auto pre = high_resolution_clock::now();

    auto duration0 = duration_cast<microseconds>(pre - start);
    cout << "Time taken for initial allocation : "<< duration0.count()/(1000.0) << " ms" << endl;

    readData(m, d, k1, blocks1, indices1, mat1, input_file1);
    readData(m, d, k2, blocks2, indices2, mat2, input_file2);
    cout <<"Reading files done!" << endl;
    auto read = high_resolution_clock::now();

    auto duration1 = duration_cast<microseconds>(read - pre);
    cout << "Time taken for reading : "<< duration1.count()/(1000.0) << " ms" << endl;

    int num_threads_per_block = 1024;
    int num_blocks = d*d/num_threads_per_block + 1;

    cout << "Number of blocks: " << num_blocks << endl;
    cout << "Going to enter the kernel" <<endl;

    //Kernel function call
    calculateResult<<<num_blocks, num_threads_per_block>>>(n, m, k1, k2, blocks1, indices1, mat1, blocks2, indices2, mat2, blocksOutput, matOutput);

    cudaError_t err = cudaGetLastError();

    if ( err != cudaSuccess )
    {
    printf("CUDA Error: %s\n", cudaGetErrorString(err));       

    // Possibly: exit(-1) if program cannot continue....
    }

    cudaDeviceSynchronize();

    auto process = high_resolution_clock::now();

    auto duration2 = duration_cast<microseconds>(process - read);
    cout << "Time taken for processing : "<< duration2.count()/(1000.0) << " ms" << endl;

    cout<<"Kernel job done !"<<endl;

    outputResult(n, m, blocksOutput, matOutput, output_file);

    auto end = high_resolution_clock::now();

    auto duration3 = duration_cast<microseconds>(end - process);
    cout << "Time taken for outputting : "<< duration3.count()/(1000.0) << " ms" << endl;

    auto duration4 = duration_cast<microseconds>(end - start);
    cout << "Total Time taken : "<< duration4.count()/(1000.0) << " ms" << endl;

    return 0;
}
