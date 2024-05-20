#include <iostream> 
#include <queue> 
using namespace std; 
#define N 6  
bool bfs(int adj[N][N],int s,int t,int parent[]) 
{ 
    char colour[N]; 
    for(int v=0;v<N;v++) 
    { 
        colour[v]='w'; 
        parent[v]=-1; 
    } 
    colour[s]='g'; 
    queue<int>q; 
    q.push(s); 
    while(!q.empty()) 
    { 
        int u=q.front(); 
        q.pop(); 
        for(int v=0;v<N;v++) 
        { 
            if(adj[u][v]>0 && colour[v]=='w') 
            { 
                parent[v]=u; 
                colour[v]='g'; 
                if(v==t) 
                    return true; 
                else 
                    q.push(v); 
            } 
        } 
        colour[u]='b'; 
 
    } 
    return false; 
} 
int ford_fulkerson(int adj[N][N],int s,int t) 
{ 
    int residual[N][N]; 
    for(int u=0;u<N;u++) 
        for(int v=0;v<N;v++) 
            residual[u][v]=adj[u][v]; 
    int max_flow=0; 
    int parent[N]; 
    while(bfs(residual,s,t,parent)) 
    { 
       int min_capacity= 9999999; 
       int v=t; 
       while(v!=s) 
       { 
            int u=parent[v]; 
            if(min_capacity>residual[u][v]) 
                min_capacity=residual[u][v]; 
            v=parent[v]; 
       } 
       v=t; 
       while (v!=s) 
       { 
            int u=parent[v]; 
            residual[u][v]-=min_capacity; 
            residual[v][u]+=min_capacity; 
            v=parent[v]; 
       } 
       max_flow+=min_capacity; 
    } 
    return max_flow;   
} 
int main() 
{ 
    int s,t; 
    int adj[N][N] 
        = { { 0, 16, 13, 0, 0, 0 }, { 0, 0, 10, 12, 0, 0 }, 
            { 0, 4, 0, 0, 14, 0 },  { 0, 0, 9, 0, 0, 20 }, 
            { 0, 0, 0, 7, 0, 4 },   { 0, 0, 0, 0, 0, 0 } }; 
    cout<<"Enter Source vertex: "; 
    cin>>s; 
    cout<<"Enter sink vertex: "; 
    cin>>t; 
    cout<<"Max Flow:"<<ford_fulkerson(adj,s,t); 
    return 0; 
} 