#ifndef NODE_H
#define NODE_H
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

namespace weakarg
{

/**
    @brief A node of the tree
*/
class Node
{
protected:
    Node * left;///<Left daughter node
    Node * right;///<Right daughter node
    Node * father;///<Father node
    double age;///<Age of the node
    int id;///<Index of the node
public:
    Node();
    Node(std::string * newick,Node * father,std::vector<Node*> * nodes);
    ~Node();
    Node* cloneSubtree() const;
    std::string newick(int p=6) const;///<Convert the node to a Newick string (precision options)
    std::string newickNoInternalLabels(int p=6) const;///<Convert the node to a Newick string (precision options)
    inline int getId() const
    {
        return id;
    }///<Returns the index of the node
    inline void setId(int i)
    {
        id=i;
    }///<Sets the index of the node
    inline double getDist() const
    {
        if (father==NULL)
            return 0;
        return father->age-age;
    }///<Returns the distance of the node to its father
    inline double getAge() const
    {
        return age;
    }///<Returns the age of the node
    inline void setAge(double a)
    {
        age=a;
    }///<Sets the age of the node
    inline Node * getLeft() const
    {
        return left;
    }///<Returns the left daughter node
    inline Node * getRight() const
    {
        return right;
    }///<Returns the right daughter node
    inline Node * getFather() const
    {
        return father;
    }///<Returns the father node
    inline void setFather(Node * f)
    {
        father=f;
    }///<Sets the father of the Node
    inline void setLeft(Node * l)
    {
        left=l;
    }///<Sets the left child of the Node
    inline void setRight(Node * r)
    {
        right=r;
    }///<Sets the right child of the Node
    inline void changeAge(double a)
    {
        age=age+a;
    }///<Changes the age of the node by adding an amount a
};

} // end namespace weakarg
#endif
