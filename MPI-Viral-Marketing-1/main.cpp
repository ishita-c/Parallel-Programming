#include <mpi.h>
#include <bits/stdc++.h>

using namespace std;

unordered_map<int, unordered_map<int, int>> adj; //adjacency list of my chunk
vector<int> degrees; // Degree of all nodes
vector<int> offsets; // offsets for input file
vector<int> owner_rank;
map<pair<int, int>, int> support;
map<pair<int, int>, vector<int>> triangles;
map<pair<int, int>, int> truss_number_for_edge;
map<set<int>, int> truss_number_for_triangle;
map<pair<int, int>, int> g;
map<pair<int, int>, unordered_map<int,int>> h;
int kmin = INT_MAX, kmax = 0;


int readFourBytes(fstream& fs){
    unsigned char a[4];
    fs.read((char *)&a, sizeof(a));
    int x = (int)a[0] | (int)a[1]<<8 | (int)a[2]<<16 | (int)a[3]<<24; 
    return x;
}

void read_graph(int num_nodes, int size, int rank, string graph_file){
    int step = num_nodes/size;
    int start_node = (num_nodes/size)*rank + min(num_nodes % size, rank);
    int end_node = (num_nodes/size)*(rank + 1) + min(num_nodes % size, rank + 1);
    
    fstream fs(graph_file, ios::in | ios::binary);
    fs.seekg(offsets[start_node], fs.beg); // start reading from your own chunk only

    for(int i=start_node;i<end_node;i++){
        int node = readFourBytes(fs);
        int deg = readFourBytes(fs);
        for (int j=0;j<deg;j++){
            int v = readFourBytes(fs);
            if(degrees[node] < degrees[v]){
                adj[node][v] = 1;
            }
            if(degrees[node] == degrees[v] and node < v){
                adj[node][v] = 1;
            }
        }
    }
    fs.close();

    // cout << "Edges in rank " << rank << ": " << count << endl;
}

void triangle_enumeration(int n, int size, int rank){
    int start_node = (n/size)*rank + min(n % size, rank);
    int end_node = (n/size)*(rank + 1) + min(n % size, rank + 1);

    int step = n/size;

    for(int i_=0;i_< (step + 1);i_++){
        vector<vector<int>> send_data(size, vector<int>());
        int u = i_ + start_node;
        if(u<end_node){
            for(auto it = adj[u].begin(); it!=adj[u].end() ;it++){
                for(auto it2 = adj[u].begin(); it2!=adj[u].end() ;it2++){
                    int v = it->first;
                    int w = it2->first;
                    if(v == w) continue;
                    if(degrees[v] > degrees[w] || ((degrees[v] == degrees[w]) and v > w)){
                        continue;
                    }
                    send_data[owner_rank[v]].push_back(u);
                    send_data[owner_rank[v]].push_back(v);
                    send_data[owner_rank[v]].push_back(w);
                }
            }
        }

        int send_counts[size];
        int sdispls[size];
        
        for(int i=0;i<size;i++) {
            send_counts[i] = send_data[i].size();
            if(i==0) {
                sdispls[i] = 0;
            } else {
                sdispls[i] = sdispls[i-1] + send_counts[i-1];
            }
            // cout << "Sending " << send_counts[i] << " from rank " << rank << " to rank " << i << endl;
        }

        int recv_counts[size];
        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

        int rdispls[size];
        for(int i=0;i<size;i++) {
            if(i==0) {
                rdispls[i] = 0;
            } else {
                rdispls[i] = rdispls[i-1] + recv_counts[i-1];
            }
        }
        int send_total = sdispls[size - 1] + send_counts[size - 1];
        int recv_total = rdispls[size - 1] + recv_counts[size - 1];

        int* send_buf = new int[send_total];
        int* recv_buf = new int[recv_total];

        int cnt = 0;
        for(int i=0;i<size;i++){
            for(int j = 0;j<send_data[i].size();j++){
                send_buf[cnt++] = send_data[i][j];
            }
        }

        MPI_Alltoallv(send_buf, send_counts, sdispls, MPI_INT,
                recv_buf, recv_counts, rdispls, MPI_INT, MPI_COMM_WORLD);

        vector<vector<int>> send_data2(size, vector<int>());

        for(int i=0;i<recv_total;i+=3){
            int u = recv_buf[i];
            int v = recv_buf[i+1];
            int w = recv_buf[i+2];

            if(adj[v].find(w)!=adj[v].end()){
                support[{v,w}]++;
                triangles[{v,w}].push_back(u);

                truss_number_for_triangle[{u,v,w}] = INT_MAX;

                send_data2[owner_rank[u]].push_back(u);
                send_data2[owner_rank[u]].push_back(v);
                send_data2[owner_rank[u]].push_back(w);
            }
        }

        memset(send_counts, 0, size * sizeof(int));
        memset(recv_counts, 0, size * sizeof(int));
        memset(sdispls, 0, size * sizeof(int));
        memset(rdispls, 0, size * sizeof(int));

        for(int i=0;i<size;i++) {
            send_counts[i] = send_data2[i].size();
            if(i==0) {
                sdispls[i] = 0;
            } else {
                sdispls[i] = sdispls[i-1] + send_counts[i-1];
            }
            // cout << "Sending " << send_counts[i] << " from rank " << rank << " to rank " << i << endl;
        }

        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

        for(int i=0;i<size;i++) {
            if(i==0) {
                rdispls[i] = 0;
            } else {
                rdispls[i] = rdispls[i-1] + recv_counts[i-1];
            }
        }
        int send_total2 = sdispls[size - 1] + send_counts[size - 1];
        int recv_total2 = rdispls[size - 1] + recv_counts[size - 1];

        delete[] send_buf;
        delete[] recv_buf;

        int* send_buf2 = new int[send_total2];
        int* recv_buf2 = new int[recv_total2];

        cnt = 0;
        for(int i=0;i<size;i++){
            for(int j = 0;j<send_data2[i].size();j++){
                send_buf2[cnt++] = send_data2[i][j];
            }
        }

        MPI_Alltoallv(send_buf2, send_counts, sdispls, MPI_INT,
                recv_buf2, recv_counts, rdispls, MPI_INT, MPI_COMM_WORLD);
        
        for(int i=0;i<recv_total2;i+=3){
            int u = recv_buf2[i];
            int v = recv_buf2[i+1];
            int w = recv_buf2[i+2];

            support[{u,v}]++;
            triangles[{u,v}].push_back(w);

            support[{u,w}]++;
            triangles[{u,w}].push_back(v);

            truss_number_for_triangle[{u,v,w}] = INT_MAX;
        }

        delete[] send_buf2;
        delete[] recv_buf2;
    }


    for(auto it = support.begin(); it!=support.end();it++){
        truss_number_for_edge[it->first] = it->second + 2;
        g[it->first] = it->second;
        kmin = min(kmin, it->second);
        kmax = max(kmax, it->second);
        for(int j=0; j<truss_number_for_edge[it->first]; j++){
            h[it->first][j] = 0;
        }
    }

    // cout << rank << ": " << support.size() << endl;
}

void truss_computation(int size, int rank){
    vector<pair<int,int>> curr;
    unordered_map<int, vector<pair<int,int>>> helper;

    for(auto it = support.begin(); it!=support.end();it++){
        helper[it->second].push_back(it->first);
    }
    int currk = kmin;
    curr = helper[kmin];
    int itn = -1;
    while(true){
        itn++;
        // cout <<"Rank "<< rank << " : " << itn << " curr=" << curr.size() << " next=";
        vector<pair<int,int>> next;
        vector<vector<int>> send_data(size, vector<int>());

        map<vector<int>, int> vis_tri;
        map<pair<int,int>, int> vis_edge; 

        map<set<int>, int> val_old_tri = truss_number_for_triangle;

        for(int i=0;i<curr.size();i++){
            int u = curr[i].first;
            int v = curr[i].second;
            if(vis_edge[{u,v}]){
                    continue;
            }else {
                vis_edge[{u,v}] = 1;
            }
            for(int j=0;j<triangles[curr[i]].size();j++){
                int w = triangles[curr[i]][j];
                if(vis_tri[{u,v,w}]){
                    continue;
                }else {
                    vis_tri[{u,v,w}] = 1;
                }

                if(adj[u].find(w)!=adj[u].end()){
                    if(truss_number_for_edge[{u,w}] < truss_number_for_edge[{u,v}]){
                        swap(v,w);
                    }
                }

            
                if(truss_number_for_edge[{u,v}] < truss_number_for_triangle[{u,v,w}]){
                    // int val_old = truss_number_for_triangle[{u,v,w}];
                    truss_number_for_triangle[{u,v,w}] = truss_number_for_edge[{u,v}];
                    int val_new = truss_number_for_triangle[{u,v,w}];

                    if((degrees[v] < degrees[w]) || ((degrees[v] == degrees[w]) and v < w)){
                        send_data[owner_rank[v]].push_back(v);
                        send_data[owner_rank[v]].push_back(w);
                        // send_data[owner_rank[v]].push_back(val_old);
                        send_data[owner_rank[v]].push_back(val_new);
                        send_data[owner_rank[v]].push_back(u);
                    } else {
                        send_data[owner_rank[w]].push_back(w);
                        send_data[owner_rank[w]].push_back(v);
                        // send_data[owner_rank[w]].push_back(val_old);
                        send_data[owner_rank[w]].push_back(val_new);
                        send_data[owner_rank[w]].push_back(u);
                    }

                    if((degrees[u] < degrees[w]) || ((degrees[u] == degrees[w]) and u < w)){
                        send_data[owner_rank[u]].push_back(u);
                        send_data[owner_rank[u]].push_back(w);
                        // send_data[owner_rank[u]].push_back(val_old);
                        send_data[owner_rank[u]].push_back(val_new);
                        send_data[owner_rank[u]].push_back(v);
                    } else {
                        send_data[owner_rank[w]].push_back(w);
                        send_data[owner_rank[w]].push_back(u);
                        // send_data[owner_rank[w]].push_back(val_old);
                        send_data[owner_rank[w]].push_back(val_new);
                        send_data[owner_rank[w]].push_back(v);
                    }
                }
            }
        }
        int send_counts[size];
        int sdispls[size];
        
        for(int i=0;i<size;i++) {
            send_counts[i] = send_data[i].size();
            if(i==0) {
                sdispls[i] = 0;
            } else {
                sdispls[i] = sdispls[i-1] + send_counts[i-1];
            }
            // cout << "Sending " << send_counts[i] << " from rank " << rank << " to rank " << i << endl;
        }

        int recv_counts[size];
        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
        int rdispls[size];
        for(int i=0;i<size;i++) {
            if(i==0) {
                rdispls[i] = 0;
            } else {
                rdispls[i] = rdispls[i-1] + recv_counts[i-1];
            }
        }
        int send_total = sdispls[size - 1] + send_counts[size - 1];
        int recv_total = rdispls[size - 1] + recv_counts[size - 1];

        int* send_buf = new int[send_total];
        int* recv_buf = new int[recv_total];

        int cnt = 0;
        for(int i=0;i<size;i++){
            for(int j = 0;j<send_data[i].size();j++){
                send_buf[cnt++] = send_data[i][j];
            }
        }

        MPI_Alltoallv(send_buf, send_counts, sdispls, MPI_INT,
                recv_buf, recv_counts, rdispls, MPI_INT, MPI_COMM_WORLD);

        map<set<int>, int> val_new_tri;

        for(int i=0; i<recv_total; i+=4){
            int u = recv_buf[i];
            int v = recv_buf[i+1];
            int w = recv_buf[i+3];
            
            if(val_new_tri.find({u,v,w})!=val_new_tri.end()){
                val_new_tri[{u,v,w}] = min(recv_buf[i+2], val_new_tri[{u,v,w}]);
            }else {
                val_new_tri[{u,v,w}] = recv_buf[i+2];
            }
        }

        map<vector<int>, int> vis;

        for(int i=0; i<recv_total; i+=4){
            int u = recv_buf[i];
            int v = recv_buf[i+1];
            int w = recv_buf[i+3];
            int val_old = val_old_tri[{u,v,w}];
            int val_new = val_new_tri[{u,v,w}];

            if(vis[{u,v,w}]){
                continue;
            }else {
                vis[{u,v,w}] = 1;
            }

            if (val_old >= truss_number_for_edge[{u,v}] and val_new < truss_number_for_edge[{u,v}])
            {
                g[{u,v}]--;
                h[{u,v}][val_new]++;
                truss_number_for_triangle[{u,v,w}] = min(truss_number_for_triangle[{u,v,w}], val_new);
            } 
            else if (val_old < truss_number_for_edge[{u,v}] and val_new < truss_number_for_edge[{u,v}])
            {
                h[{u,v}][val_old]--;
                h[{u,v}][val_new]++;
                truss_number_for_triangle[{u,v,w}] = min(truss_number_for_triangle[{u,v,w}], val_new);
            }
            if(g[{u,v}] < truss_number_for_edge[{u,v}] - 2){
                truss_number_for_edge[{u,v}]--;
                g[{u,v}] = g[{u,v}] + h[{u,v}][truss_number_for_edge[{u,v}]];
                next.push_back({u,v});
            }
        }

        delete[] send_buf;
        delete[] recv_buf;

        currk += 1;
        for(int j=0;j<helper[currk].size();j++){
            next.push_back(helper[currk][j]);
        }

        // cout << next.size() << endl;

        int all_sum;
        int my_num = next.size();
        MPI_Allreduce(&my_num, &all_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(all_sum == 0) break;
        else curr = next;
    }
}

void dfs(int i, int k, vector<int> &vis, vector<int> &component, vector<vector<pair<int,int>>> &adj_dfs){
    vis[i] = 1;
    component.push_back(i);
    for(auto p: adj_dfs[i]){
        if(!vis[p.first] and p.second >= k + 2){
            dfs(p.first, k, vis, component, adj_dfs);
        }
    }
}

void output_result(int n, int size, int rank, int startk, int endk, string output_path, int verbose){

    if (verbose == 0){
        int mxks[size];
        int mxk = 0;
        int root = 0;
        for(auto it=truss_number_for_edge.begin();it!=truss_number_for_edge.end();it++){
            mxk = max(mxk, it->second);
        }
        MPI_Gather(&mxk, 1, MPI_INT, mxks, 1, MPI_INT, root, MPI_COMM_WORLD);
        if (rank == 0){
            for(int i=0;i<size;i++){
                mxk = max(mxk, mxks[i]);
            }
            fstream fs(output_path, ios::out);
            for(int k=startk;k<=endk;k++){
                fs << (k <= mxk - 2) << " ";
            }
            fs.close();
        }

    }
    else {
        int counts[size];
        int sendcount = truss_number_for_edge.size()*3;
        int root = 0;
        int displs[size];

        MPI_Gather(&sendcount, 1, MPI_INT, counts, 1, MPI_INT, root, MPI_COMM_WORLD);

        if (rank == 0) {
            displs[0] = 0;
            for (int i = 1; i < size; i++) {
                displs[i] = displs[i - 1] + counts[i - 1];
            }
        }

        int recv_total = displs[size-1] + counts[size-1];

        if(rank !=0 ) recv_total = 0;

        int* sendbuf = new int[sendcount];
        int* recvbuf = new int[recv_total];

        int cnt = 0;
        for(auto it=truss_number_for_edge.begin();it!=truss_number_for_edge.end();it++){
            sendbuf[cnt++] = it->first.first;
            sendbuf[cnt++] = it->first.second;
            sendbuf[cnt++] = it->second;
        }

        MPI_Gatherv(sendbuf, sendcount, MPI_INT, recvbuf, counts, displs, MPI_INT, root, MPI_COMM_WORLD);

        if (rank == 0) {
            int mxk = 0;
            vector<vector<pair<int,int>>> adj_dfs(n);
            for(int i=0;i<recv_total;i+=3){
                adj_dfs[recvbuf[i]].push_back({recvbuf[i+1], recvbuf[i+2]});
                adj_dfs[recvbuf[i+1]].push_back({recvbuf[i], recvbuf[i+2]});
                mxk = max(mxk, recvbuf[i+2]);
            }

            fstream fs(output_path, ios::out);

            for(int k=startk;k<=endk;k++){
                fs << (k <= mxk - 2) << endl;
                if(k > mxk - 2) continue;
                vector<vector<int>> components;
                vector<int> vis(n,0);
                for(int i=0;i<n;i++){
                    if(!vis[i]){
                        vector<int> component;
                        dfs(i, k, vis, component, adj_dfs);
                        if(component.size() > k+1){
                            components.push_back(component);
                        }
                    }
                }
                fs << components.size() << endl;
                for(int i=0;i<components.size();i++){
                    // sort(components[i].begin(), components[i].end()); // Not necessary, can be commented for performance
                    for(int j=0;j<components[i].size();j++){
                        fs << components[i][j];
                        if(j != components[i].size() - 1) fs << " ";
                    }
                    fs << endl;
                }
            }
            fs.close();
        }
        delete[] sendbuf;
        delete[] recvbuf;
    }
}

int main(int argc, char* argv[]){
    cout << endl;

    // Default values for the command line arguments
    int taskid = 1;
    string input_path = "";
    string header_path = "";
    string output_path = "";
    int verbose = 0;
    int startk = -1;
    int endk = -1;

    // Loop through the command line arguments and parse them
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg.find("--taskid=") == 0) {
            taskid = std::stoi(arg.substr(9));
        } else if (arg.find("--inputpath=") != string::npos) {
            input_path = arg.substr(12);
        } else if (arg.find("--headerpath=") != string::npos) {
            header_path = arg.substr(13);
        } else if (arg.find("--outputpath=") != string::npos) {
            output_path = arg.substr(13);
        } else if (arg.find("--verbose=") != string::npos) {
            verbose = stoi(arg.substr(10));
        } else if (arg.find("--startk=") != string::npos) {
            startk = stoi(arg.substr(9));
        } else if (arg.find("--endk=") != string::npos) {
            endk = stoi(arg.substr(7));
        }
    }

    if(input_path == "" || header_path == "" || output_path == "" || startk == -1 || endk == -1){
        cout << "Please correct your arguements !" << endl;
    }

    // cout << "taskid = " << taskid << endl;
    // cout << "input_path = " << input_path << endl;
    // cout << "header_path = " << header_path << endl;
    // cout << "output_path = " << output_path << endl;
    // cout << "verbose = " << verbose << endl;
    // cout << "startk = " << startk << endl;
    // cout << "endk = " << endk << endl;

    fstream fs_input(input_path,std::ios::in | std::ios::binary);
    int n = readFourBytes(fs_input);
    int m = readFourBytes(fs_input);

    fstream fs_header(header_path,std::ios::in | std::ios::binary);
    offsets.resize(n);
    for(int i=0;i<n;i++){
        offsets[i] = readFourBytes(fs_header);
    }
    fs_header.close();

    degrees.resize(n);
    for(int i=0;i<n;i++){
        fs_input.seekg(offsets[i] + 4, fs_input.beg);
        degrees[i] = readFourBytes(fs_input);
    }
    fs_input.close();

    int rank, size;

    //Starting MPI pipeline
    MPI_Init(NULL, NULL);
    
    // Extracting Rank and Processor Count
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Assign owners
    owner_rank.resize(n);
    for(int i=0;i<size;i++){
        int start_node = (n/size)*i + min(n % size, i);
        int end_node = (n/size)*(i + 1) + min(n % size, i + 1);
        for(int j=start_node;j<end_node;j++){
            owner_rank[j] = i;
        }
        if(rank == 0){
            cout << "Range for rank " << i << ": [" << start_node << ", " <<  end_node << ")" << endl; 
        }
    }

    // Read Graph
    read_graph(n, size, rank, input_path);

    cout << "Graph reading done!" << endl;


    // Preprocessing - Triangle Enumeration
    triangle_enumeration(n, size, rank);

    // for(auto it = support.begin();it!=support.end();it++){
    //     cout << rank << " : " << it->first.first << "," << it->first.second << " = " << it->second << endl;
    // }

    cout << "Triangle Enumeration done!" << endl;


    // Truss Computation
    truss_computation(size, rank);

    // for(auto it = truss_number_for_edge.begin();it!=truss_number_for_edge.end();it++){
    //     cout << rank << " : " << it->first.first << "," << it->first.second << " = " << it->second << endl;
    // }

    cout << "Truss Computation done!" << endl;

    // Output result
    output_result(n, size, rank, startk, endk, output_path, verbose);

    cout << "Outputting Result done!" << endl;

    MPI_Finalize();
    
    return 0;
}