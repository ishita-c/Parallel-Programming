#include <bits/stdc++.h>
#include <omp.h>
#include "library.hpp"

using namespace std;

int MAX_VAL = pow(2,16) - 1;

int readFourBytes(fstream& fs){
    unsigned char a[4];
    fs.read((char *)&a, sizeof(a));
    int x = (int)a[0] | (int)a[1]<<8 | (int)a[2]<<16 | (int)a[3]<<24; 
    return x;
}

int readOneByte(fstream& fs){
    unsigned char a[1];
    fs.read((char *)&a, sizeof(a));
    int x = (int)a[0];
    return x;
}

void writeFourBytes(int x, fstream& fs){
    char* buffer;
    buffer = (char*)(&x);

    for(int i = 0; i < 4; i++)
        fs.write(buffer + i, sizeof(char));
}

void writeTwoBytes(int x, fstream& fs){
    char* buffer;
    buffer = (char*)(&x);

    for(int i = 0; i < 2; i++)
        fs.write(buffer + i, sizeof(char));
}

void readData(int m, int d, vector<pair<int, vector<int>>> &data, string input_file){
    int k = data.size();
    int numt = omp_get_num_threads();
    int blocksize = m*m + 8;
    int step = k/numt;

    for(int i_ = 0 ; i_ < numt ; i_++){
        #pragma omp task shared(data)
        {
            int start = i_*step;
            int end = (i_==(numt-1)?k:(i_+1)*step);
            fstream fs(input_file, ios::in | ios::binary);
            fs.seekg(blocksize*start + 12, std::ios::beg);
            for(int j_ = start ; j_ < end ; j_++){
                int i = readFourBytes(fs);
                int j = readFourBytes(fs);
                vector<int> temp(m*m);
                for(int k_ = 0; k_ < m*m; k_++){
                    int x = readOneByte(fs);
                    temp[k_] = x;
                }
                data[j_] = {i*d + j, temp};
            }
            fs.close();
        }
    }
    #pragma omp taskwait
    cout << "Reading Data done!"<<endl;
}

void completeData(int n, int m, vector<pair<int, vector<int>>> &data, vector<int>& stats, unordered_map<int, vector<pair<int, vector<int>>>>& pre){
    int k = data.size();
    int d = n/m;
    for(int i = 0; i < k; i++){
        int x = data[i].first/d;
        int y = data[i].first%d;
        stats[x+1]++;
        if(x == y) {
            for(int p = 0; p<m; p++){
                for(int q = 0; q<p; q++){
                    data[i].second[p*m + q] = data[i].second[q*m + p];
                }
            }
            pre[x].push_back(data[i]);
            continue;
        }
        pre[x].push_back(data[i]);
        stats[y+1]++;
        vector<int> temp(m*m);
        for(int j=0;j<m*m;j++){
            temp[(j%m)*m + j/m] = data[i].second[j];
        }
        data.push_back({y*d + x, temp});
        pre[y].push_back({y*d + x, temp});
    }
    for(int i=1;i<=d;i++){
        stats[i] += stats[i-1];
    }
    sort(data.begin(), data.end());
    cout << "Preprocessing Data done!"<<endl;
}

void writeResult(int m, int d, vector<unordered_map<int, vector<int>>> &res, string output_file){
    int numt = omp_get_num_threads();
    int blocksize = 2*m*m + 8;
    // int step = k/numt;
    int n = d*m;

    vector<int> ranges(numt + 1, 0);
    for(int i=0;i<numt;i++) {
        ranges[i + 1] = res[i].size();
        if(i>0) ranges[i + 1] += ranges[i];
    }

    int k = ranges[numt];

    fstream fs(output_file, ios::out | ios::binary);
    writeFourBytes(n, fs);
    writeFourBytes(m, fs);
    writeFourBytes(k, fs);
    fs.close();

    for(int i_ = 0 ; i_ < numt ; i_++){
        #pragma omp task shared(res)
        {
            int start = ranges[i_];
            int end = ranges[i_ + 1];
            unordered_map <int, vector<int>>::iterator it = res[i_].begin();

            fstream fs(output_file, ios::in | ios::binary | ios::out);
            fs.seekp(blocksize*start + 12, ios::beg);

            for(auto j=start; j<end; j++, it++){
                int idx = it->first;
                writeFourBytes(idx/d, fs);
                writeFourBytes(idx%d, fs);
                for(int k_ = 0; k_ < m*m; k_++){
                    writeTwoBytes(it->second[k_], fs);
                }
            }
            fs.close();
        }
    }
    #pragma omp taskwait

    fs.close();
    cout<< "Writing Result done!"<<endl;
    cout << "Output k = " << k << endl;
}

bool innerBlocks(int m, vector<int>& a, vector<int>& b, vector<int>& c){
    bool flag = false;
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            int val = 0;
            for(int k=0;k<m;k++){
                val = Outer(val, Inner(a[i*m + k], b[k*m + j]));
            }
            if(val > 0){
                flag = true;
                c[i*m + j] = min(val, MAX_VAL);
            } 
        }
    }
    return flag;
}

void outerBlocks(int m, vector<int>& a, vector<int>& b){
    for(int i=0;i<m*m;i++){
        a[i] = min(Outer(a[i], b[i]), MAX_VAL);
    }
}

void calculateResultEfficient(int m, int d, vector<pair<int, vector<int>>> &data, vector<int> &stats, unordered_map<int, vector<pair<int, vector<int>>>>& pre, vector<unordered_map<int, vector<int>>>& res){
    int numt = omp_get_num_threads();
    int step = d/numt;

    for(int i_ = 0; i_ < numt; i_++){
        #pragma omp task shared(res)
        {
            int start = stats[i_*step];
            int end = stats[(i_==(numt-1)?d:(i_+1)*step)];
            unordered_map<int, int> mp;
            vector<int> c(m*m, 0);

            for(int i=start;i<end;i++){
                int x1 = data[i].first/d;
                int y1 = data[i].first%d;
                for(int j=0;j< int(pre[y1].size());j++){
                    // int x2 = pre[y1][j].first/d;
                    int y2 = pre[y1][j].first%d;

                    if(x1 > y2) continue;

                    c.assign(m*m, 0);
                    if(innerBlocks(m, data[i].second, pre[y1][j].second, c)){
                        if(mp[x1*d + y2] == 1){
                            outerBlocks(m, res[i_][x1*d + y2], c);
                        }else {
                            res[i_][x1*d + y2] = c;
                            mp[x1*d + y2] = 1;
                        }
                    }
                }
            }
        }
    }
    #pragma omp taskwait

    cout<< "Calculating Result done! " << endl;
}


int main(int argc, char* argv[]){

    string input_file = argv[1];
    string output_file = argv[2];
    
    fstream fs(input_file, ios::in | ios::binary);

    int n = readFourBytes(fs);
    int m = readFourBytes(fs);
    int k = readFourBytes(fs);
    fs.close();

    // dimension of block matrix
    int d = n/m;

    cout << "n = " << n << endl;
    cout << "m = " << m << endl;
    cout << "k = " << k << endl;

    int n_threads = 8;  

    if (k > pow(2,13)) n_threads = 16;    

    vector<pair<int, vector<int>>> data(k);
    vector<int> stats(d+1,0);
    unordered_map<int, vector<pair<int, vector<int>>>> pre;
    vector<unordered_map<int, vector<int>>> res(n_threads);

    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp single
        { 
            auto start = chrono::high_resolution_clock::now();

            readData(m, d, data, input_file);
            auto read = chrono::high_resolution_clock::now();

            completeData(n, m, data, stats, pre);
            auto complete = chrono::high_resolution_clock::now();

            calculateResultEfficient(m, d, data, stats, pre, res);
            auto calculate = chrono::high_resolution_clock::now();

            writeResult(m, d, res, output_file);
            auto write = chrono::high_resolution_clock::now();

            auto duration = chrono::duration_cast<chrono::microseconds>(read - start);
            cout << "Read time: "<<duration.count()/1000 << endl;

            auto duration2 = chrono::duration_cast<chrono::microseconds>(complete - read);
            cout << "Preprocessing time: "<<duration2.count()/1000 << endl;

            auto duration3 = chrono::duration_cast<chrono::microseconds>(calculate - complete);
            cout << "Calculate time: "<<duration3.count()/1000 << endl;

            auto duration4 = chrono::duration_cast<chrono::microseconds>(write - calculate);
            cout << "Write time: "<<duration4.count()/1000 << endl;

            auto duration5 = chrono::duration_cast<chrono::microseconds>(write - start);
            cout << "Total time: "<<duration5.count()/1000 << endl;
        } 
    }
    return 0;
}