#include "node.h"
#include "slotallocator.h"
#include <sstream>


using namespace std;
namespace weakarg
{

Node::Node()
{
    left=NULL;
    right=NULL;
    father=NULL;
    age=0.0;
    id=0;
}

Node::Node(string * newick,Node * father,vector<Node*> * nodes)
{
	age=0.0;
	id=0;
    this->father=father;
    //Find central comma if it exists
    int depth=0;
    int found=0;
    for (unsigned int i=0;i<newick->length();i++)
    {
        if (newick->at(i)=='(')
        {
            depth++;
            continue;
        }
        if (newick->at(i)==')')
        {
            depth--;
            continue;
        }
        if (newick->at(i)==',' && depth==1)
        {
            found=i;
            break;
        }
    }
    //Find last ':' and read dist
    int found2=0;
    depth=0;
    for (unsigned int i=0;i<newick->length();i++)
    {
        if (newick->at(i)=='(')
            depth++;
        if (newick->at(i)==':')
            found2=i;
        if (newick->at(i)==')')
        {
            depth--;
            if (depth==0)
                found2=0;
        }
    }
    if (found2>0)
    {
        istringstream input(newick->substr(found2+1,newick->length()));
        input >> age;
    }
    else
    {
        age=0.0;
        found2=newick->length();
    }

    if (found==0)
    {
        //Leaf
        left=NULL;
        right=NULL;
        istringstream input(*newick);
        input >> id;
        nodes->push_back(this);
    }
    else
    {
        //Internal node
        id=-1;
        string leftStr =newick->substr(1,found-1);
        string rightStr=newick->substr(found+1,found2-2-found);
        left=SlotAllocator<Node>::GetSlotAllocator().Allocate();	// allocate storage
        right=SlotAllocator<Node>::GetSlotAllocator().Allocate();
        new (left) Node(&leftStr,this,nodes);	// use in-place constructor
        new (right) Node(&rightStr,this,nodes);
        nodes->push_back(this);
    }
}

Node* Node::cloneSubtree() const
{
	Node* newnode = SlotAllocator<Node>::GetSlotAllocator().Allocate();
	newnode->age = age;
	newnode->id = id;
	newnode->left = NULL;
	newnode->right = NULL;
	if(left!=NULL)
	{
		newnode->left = left->cloneSubtree();
		newnode->left->father = newnode;
	}
	if(right!=NULL)
	{
		newnode->right = right->cloneSubtree();
		newnode->right->father = newnode;
	}
	newnode->father = father;	// ensures that root node gets null father
	return newnode;
}

string Node::newick(int p) const
{
    ostringstream idstream;
    ostringstream diststream;
    diststream.setf(ios::scientific, ios::floatfield);
    diststream.precision(p);
    idstream<<id;
    diststream<<getDist();
    if (left==NULL)
        return idstream.str()+":"+diststream.str();
    else return "("+left->newick(p)+","+right->newick(p)+")"+idstream.str()+":"+diststream.str();
}

string Node::newickNoInternalLabels(int p) const
{
    ostringstream idstream;
    ostringstream diststream;
    diststream.setf(ios::scientific, ios::floatfield);
    diststream.precision(p);
    idstream<<id;
    diststream<<fixed<<setprecision (9)<<getDist();
    if (left==NULL)
        return idstream.str()+":"+diststream.str();
    else         return "("+left->newickNoInternalLabels(p)+","+right->newickNoInternalLabels(p)+")"+":"+diststream.str();
}

Node::~Node()
{
    if (left )
    	SlotAllocator<Node>::GetSlotAllocator().Free(left );
    if (right)
    	SlotAllocator<Node>::GetSlotAllocator().Free(right);
    left=NULL;
    right=NULL;
}


} // end namespace weakarg
