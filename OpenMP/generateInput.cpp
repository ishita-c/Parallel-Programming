#include <bits/stdc++.h>
using namespace std;

void writeFourBytes(int x, fstream& fs){
    char* buffer;
    buffer = (char*)(&x);

    for(int i = 0; i < 4; i++)
        fs.write(buffer + i, sizeof(char));
}

void writeOneByte(int x, fstream& fs){
    char* buffer;
    buffer = (char*)(&x);
    fs.write(buffer, sizeof(char));
}

int main(int argc, char* argv[]){

    int n = 100000;
    int m = 25;
    int k = pow(2, atoi(argv[1]));
    string output_file = "test_" + string(argv[1]);
    fstream fs(output_file, ios::out | ios::binary);

    int d = n/m;

    writeFourBytes(n, fs);
    writeFourBytes(m, fs);
    writeFourBytes(k, fs);

    unordered_map<int, int > mp;

    for(int i=0;i<k;i++){
        int x,y;
        while(true){
            x = rand()%d;
            y = rand()%d;
            if(y >= x and mp.find(x*d + y) == mp.end()){
                mp[x*d + y] = 1;
                break;
            }
        }
        writeFourBytes(x, fs);
        writeFourBytes(y, fs);
        for(int j=0;j<m*m;j++){
            int val = rand()%10;
            writeOneByte(val, fs);
        }
    }

    fs.close();
    return 0;
}