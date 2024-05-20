#include <iostream> 
using namespace std; 
class BST; 
class node{ 
    int data,tag; 
    node *right,*left; 
    public: 
        node() 
        { 
            right=NULL; 
            left=NULL; 
        } 
        node(int data) 
        { 
            this->data=data; 
            right=NULL; 
            left=NULL; 
            this->tag=0; 
        } 
        friend class BST; 
}; 
class BST{ 
    node *root; 
    public : 
        node* append(node* root,int data); 
        void Linkinversion_inorder(node* pres); 
        void Linkinversion_preorder(node* pres); 
        void Linkinversion_postorder(node* pres); 
}; 
node* BST :: append(node* root,int data) 
{ 
    if(root==NULL) 
    { 
        return new node(data); 
    } 
    if(data>root->data) 
        root->right=append(root->right,data); 
    else  if(data<root->data)   
        root->left=append(root->left,data); 
    return root; 
} 
void BST::Linkinversion_inorder(node* pres) 
{ 
    cout<<endl<<"INORDER: "; 
    node* prev=NULL; 
    node* temp=NULL; 
    int flag=0; 
    if(pres==NULL) 
        return ; 
    do{ 
        if(pres->left!=NULL && flag==0) //descend left 
        {                                         
            pres->tag=0;                             
            temp=pres->left; 
            pres->left=prev; 
            prev=pres; 
            pres=temp; 
            flag=0; 
        } 
        else{ 
            if(flag!=2)  
            { 
                cout<<pres->data<<" "; 
                if(prev==NULL && pres->right==NULL) 
                    flag=2; 
            } 
            if(pres->right!=NULL && flag!=2) //descend right 
            { 
                pres->tag=1; 
                temp=pres->right ; 
                pres->right=prev; 
                prev=pres; 
                pres=temp; 
                flag=0; 
            } 
            else if(prev!=NULL && prev->tag==0) //ascend left 
            { 
                prev->tag=0; 
                temp=prev->left;  
                prev->left=pres; 
                pres=prev; 
                prev=temp; 
                flag=1; 
            } 
            else if(prev!=NULL && prev->tag==1) //ascend right 
            { 
                prev->tag=0; 
                temp=prev->right; 
                prev->right=pres; 
                pres=prev; 
                prev=temp; 
                flag=2; 
            } 
        } 
    }while(prev!=NULL || flag!=2); 
} 
void BST :: Linkinversion_preorder(node* pres) 
{ 
    cout<<endl<<"PREORDER: "; 
    node* prev=NULL; 
    node* temp=NULL; 
    int flag=0; 
    if(pres==NULL) 
        return ; 
    do{ 
        if(flag==0) 
        { 
            cout<<pres->data<<" "; 
            if(prev==NULL && pres->right==NULL && pres->left==NULL) 
                flag=2; 
 
        } 
        if(pres->left!=NULL && flag==0) //descend left 
        { 
                pres->tag=0; 
                temp=pres->left; 
                pres->left=prev; 
                prev=pres; 
                pres=temp; 
                flag=0; 
        } 
        else{ 
                if(pres->right!=NULL && flag!=2) //descend right 
                { 
                    pres->tag=1; 
                    temp=pres->right ; 
                    pres->right=prev; 
                    prev=pres; 
                    pres=temp; 
                    flag=0; 
                } 
                else if(prev!=NULL && prev->tag==0) //ascend left 
                { 
                    prev->tag=0; 
                    temp=prev->left;  
                    prev->left=pres; 
                    pres=prev; 
                    prev=temp; 
                    if(prev==NULL && pres->right==NULL) 
                        flag=2; 
                    else 
                        flag=1; 
                } 
                else if(prev!=NULL && prev->tag==1) //ascend right 
                { 
                    prev->tag=0; 
                    temp=prev->right; 
                    prev->right=pres; 
                    pres=prev; 
                    prev=temp; 
                    flag=2; 
                } 
            } 
    }while(prev!=NULL || flag!=2); 
} 
void BST :: Linkinversion_postorder(node* pres) 
{ 
    cout<<endl<<"POSTORDER: "; 
    node* prev=NULL; 
    node* temp=NULL; 
    int flag=0; 
    if(pres==NULL) 
        return ;     
    do{ 
        if(pres->left!=NULL && flag==0) //descend left 
        { 
            pres->tag=0; 
            temp=pres->left; 
            pres->left=prev; 
            prev=pres; 
            pres=temp; 
            flag=0; 
        } 
        else{ 
            if(pres->right!=NULL && flag!=2) //descend right 
            { 
                pres->tag=1; 
                temp=pres->right ; 
                pres->right=prev; 
                prev=pres; 
                pres=temp; 
                flag=0; 
            }   
            else if(prev!=NULL && prev->tag==0) //ascend left 
            { 
                prev->tag=0; 
                temp=prev->left;  
                prev->left=pres; 
                pres=prev; 
                prev=temp; 
                flag=1; 
            } 
            else if(prev!=NULL && prev->tag==1) //ascend right 
            { 
                prev->tag=0; 
                temp=prev->right; 
                prev->right=pres; 
                pres=prev; 
                prev=temp; 
                flag=2; 
            }           
        } 
        if((pres->right==NULL&&flag==1)||(pres->left==NULL&&pres->right==NULL&&flag==0)||flag==2) 
        { 
            cout<<pres->data<<" "; 
            flag=2; 
        }  
        
    }while(prev!=NULL || flag!=2); 
} 
int main() 
{ 
    BST b; 
    node* root=NULL; 
    root=b.append(root,7); 
    b.append(root,5); 
    b.append(root,9); 
    b.append(root,4); 
    b.append(root,6); 
    b.append(root,41); 
    b.append(root,449); 
    b.append(root,54); 
    b.Linkinversion_inorder(root); 
    b.Linkinversion_preorder(root); 
    b.Linkinversion_postorder(root); 
    return 0; 
}