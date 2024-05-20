#include <iostream> 
#include <cmath> 
using namespace std; 
void binaryCounter(int a[],int k) 
{ 
    int i=0; 
    while(i<k && a[i]==1) 
    { 
        a[i]=0; 
        i++; 
    } 
    if(i<k) 
        a[i]=1; 
} 
void display(int a[],int k) 
{ 
    for(int i=k-1;i>=0;i--) 
        cout<<a[i]<<" "; 
} 
int main() 
{ 
    int k; 
    cout<<"Enter no of bit :"; 
    cin>>k; 
    int *a=new int[k]; 
    for(int i=0;i<k;i++) 
        a[i]=0; 
    int n=pow(2,k)-1; 
    display(a,k); 
    for(int i=0;i<n;i++) 
    { 
        cout<<endl; 
        binaryCounter(a,k); 
        display(a,k); 
    } 
    return 0; 
} 