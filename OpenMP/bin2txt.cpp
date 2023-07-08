#include <bits/stdc++.h>

using namespace std;

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

int main(int argc, char* argv[]){

    string input_file = argv[1];

    fstream fsin(input_file, ios::in | ios::binary);

    int n = readFourBytes(fsin);
    int m = readFourBytes(fsin);
    int k = readFourBytes(fsin);

    int d = n/m;

    map<int, vector<int>> mp;

    for(int i=0;i<k;i++){
        int x = readFourBytes(fsin);
        int y = readFourBytes(fsin);
        for(int j=0;j<m*m;j++){
            int val = readTwoBytes(fsin);
            mp[x*d + y].push_back(val);
        }
    }

    fsin.close();

    fstream fsout(input_file + ".txt", ios::out);

    fsout << "n = " << n << endl;
    fsout << "m = " << m << endl;
    fsout << "k = " << k << endl;

    int cnt = 0;

    for(auto i = mp.begin(); i != mp.end(); i++){
        int x = i->first/d;
        int y = i->first%d;
        fsout << cnt++ <<") "<< x << "," << y << " : ";
        for(int j=0;j<m*m;j++){
            fsout << i->second[j] << " ";
        }
        fsout<<endl;
    }

    fsout.close();
    return 0;
}