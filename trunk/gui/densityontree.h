#ifndef DENSITYONTREE_H
#define DENSITYONTREE_H
//
#include "gelmanrubinimpl.h"
//
class DensityOnTree  
{

public:
	DensityOnTree(Tree*tree,double timeScale);
	void display(QPainter * painter,vector<double>*x,vector<double>*y,double scale,QBrush brush);
	void add(int edge ,double age,double relcon=1.0);
	void smooth();
	//inline int maximum(){int max=0;for (unsigned int i=0;i<d.size();i++) if (max<d[i]) max=d[i];return max;}
	inline void setTree(Tree*tree) {this->tree=tree;}
	void add2(RecTree *treefrom,int i,bool goingto,double relcon=1.0);///<Adds a projection onto the tree

protected:
	Tree * tree;
	vector <double> d;
	double above;
	double precision;
	double timeScale;
	
};
#endif
