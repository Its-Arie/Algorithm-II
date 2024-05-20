#include <iostream> 
using namespace std; 
class stack{ 
    private: 
        int *arr; 
        int top; 
        int capacity; 
    public: 
        void init(int size) 
        { 
            arr= new int[size]; 
            top=-1; 
            capacity=size; 
        } 
        bool isfull() 
        { 
            return top==capacity-1; 
        } 
        bool isempty() 
        { 
            return top==-1; 
        } 
        void push(int data) 
        { 
            if(!isfull()) 
                arr[++top]=data; 
            else 
                cout<<"Stack is overflow"<<endl; 
        } 
        void pop() 
        { 
            if(!isempty()) 
                cout<<arr[top--]<<" is popped."<<endl; 
            else 
                cout<<"Stack is underflow."<<endl; 
        } 
        void multipop(int n) 
        { 
            for(int i=0;i<n;i++) 
                pop(); 
        } 
}; 
int main() 
{ 
    int ch; 
    stack s; 
    s.init(6); 
    do{ 
        cout<<"-----------------------------------"<<endl; 
        cout<<"Enter 1 for push: "<<endl; 
        cout<<"Enter 2 for pop: "<<endl; 
        cout<<"Enter 3 for multipop: "<<endl; 
        cout<<"Enter 4 for exit: "<<endl; 
        cout<<"Enter choice : "; 
        cin>>ch; 
        switch (ch) 
        { 
        case 1: 
            int x,N; 
            cout<<"How many elements do you want to push: "; 
            cin>>N; 
            for(int i=0;i<N;i++) 
            { 
                cout<<"Enter number: "; 
                cin>>x; 
                s.push(x); 
            } 
            break; 
        case 2: 
            s.pop(); 
            break; 
        case 3: 
            int n; 
            cout<<"Enter how many elements do you want to delete: "; 
            cin>>n; 
            s.multipop(n); 
            break; 
        } 
         
    }while(ch!=4); 
    return 0; 
} 