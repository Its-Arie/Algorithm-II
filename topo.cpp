#include <iostream>
#include <bits/stdc++.h>
using namespace std;
class Graph;
class Stack
{
    int *s;
    int top;
    int capacity;
    bool is_full();
    bool is_empty();
    void push(int a);
    int pop();

public:
    Stack()
    {
        capacity = 100;
        s = new int[capacity];
        top = -1;
    }
    friend class Graph;
};

bool Stack::is_full()
{
    return top == capacity - 1;
}


bool Stack::is_empty()
{
    return top == -1;
}

void Stack::push(int a)
{
    if (!is_full())
    {
        s[++top] = a;
    }
}

int Stack::pop()
{
    if (!is_empty())
    {
        return s[top--];
    }
}

class Graph
{
    int vertex;
    int **arr;
    char *color;
    Stack s;
    void visit(int v);

public:
    Graph(int vertex)
    {
        this->vertex = vertex;
        arr = new int *[vertex];
        for (int i = 0; i < vertex; i++)
        {
            arr[i] = new int[vertex];
        }
        color = new char[vertex];
    }
    void input();
    void topological_sort();

    void display();
};

void Graph::display()
{
    while (!s.is_empty())
    {
        cout << s.pop() << "\t";
    }
}

void Graph::visit(int v)
{
    color[v] = 'g';
    for (int i = 0; i < vertex; i++)
    {
        if (arr[v][i] == 1 && color[i] == 'w')
        {
            visit(i);
        }
    }
    color[v] = 'b';
    s.push(v);
}

void Graph::topological_sort()
{
    for (int i = 0; i < vertex; i++)
    {
        if (color[i] == 'w')
        {
            visit(i);
        }
    }
}

void Graph::input()
{
    cout << "input" << endl;
    for (int i = 0; i < vertex; i++)
    {
        for (int j = 0; j < vertex; j++)
        {
            cin >> arr[i][j];
        }
    }

    for (int i = 0; i < vertex; i++)
    {
        color[i] = 'w';
    }
}

int main()
{
    int vertex;
    cout << "vertex : ";
    cin >> vertex;
    Graph a(vertex);
    a.input();
    a.topological_sort();
    a.display();
    return 0;
}