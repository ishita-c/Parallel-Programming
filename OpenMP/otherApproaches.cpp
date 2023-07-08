#include <bits/stdc++.h>
#include <omp.h>
#include "library.hpp"

using namespace std;

void calculateResultBrute(int m, int d, vector<pair<int, vector<int>>> &data, vector<pair<int, vector<int>>> &res){
    // n^3 algorithm
    int n = m*d;
    int z = data.size();
    map<int, vector<int>> dict;
    for(int i=0;i<z;i++){
        int idx = data[i].first;
        dict[idx] = data[i].second;
    }
    
    for(int i = 0; i < d; i++){
        for(int j = i; j < d; j++){
            bool flag = false;
            vector<int> temp(m*m);
            for(int k = 0;k <m*m;k++){
                int x = i*m + k/m;
                int y = j*m + k%m;
                int val = 0;
                for(int p = 0; p < n; p++){
                    // val += data[x][p] * data[p][y]
                    int p1_idx1, p1_idx2, p2_idx1, p2_idx2;

                    if(x <= p){
                        p1_idx1 = (x/m)*d + p/m;
                        p1_idx2 = (x%m)*m + p%m;
                    } else {
                        p1_idx1 = (p/m)*d + x/m;
                        p1_idx2 = (p%m)*m + x%m;
                    }
                    if (p <= y) {
                        p2_idx1 = (p/m)*d + y/m;
                        p2_idx2 = (p%m)*m + y%m;
                    } else {
                        p2_idx1 = (y/m)*d + p/m;
                        p2_idx2 = (y%m)*m + p%m;
                    }

                    int val1 = 0, val2 = 0;

                    if(dict.count(p1_idx1)){
                        val1 = dict[p1_idx1][p1_idx2];
                    }
                    if(dict.count(p2_idx1)){
                        val2 = dict[p2_idx1][p2_idx2];
                    }

                    val = Outer(val, Inner(val1, val2));
                    //cout << val << " " << val1 << " " << val2 << endl; 
                }
                temp[k] = val;
                if(val!=0) flag = true;
            }
            if(flag){
                res.push_back({i*d + j, temp});
            }
        }
    }
    cout<< "Calculating Result done! k = "<< res.size() << endl;
}