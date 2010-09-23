#include "densityontree.h"
//
DensityOnTree::DensityOnTree(Tree * tree,double timeScale) 
{
	this->timeScale=timeScale;
	precision=0.005*tree->getNode(tree->getN()*2-2)->getAge();
	this->tree=tree;
	int siz=0;
	for (int i=0;i<tree->getN()*2-2;i++) siz+=floor(tree->getDist(i)/precision);
	siz+=max(0.0,floor((timeScale-tree->getNode(tree->getN()*2-2)->getAge())/precision));
	d=vector<double>(siz,0.0);
	above=0.0;
}

void DensityOnTree::add(int edge,double age,double relcon)
{
int which=0;
int n=tree->getN();
for (int e=0;e<edge;e++) which+=floor(tree->getDist(e)/precision);
int end=which;
which+=floor(age/precision);
if (edge<n*2-2) 
{end+=floor(tree->getDist(edge)/precision);if (which>=end) return;}
else 
{end+=floor((timeScale-tree->getNode(n*2-2)->getAge())/precision);if (which>=end) {above++;return;};}
if(which>=d.size()) {d[d.size()-1]+=relcon;
}else d[which]+=relcon;  //*** See below.  There is a danger that something is wrong here because when mapping to a consensus tree sometimes end>d.size....
}

void DensityOnTree::smooth()
{
	int n=tree->getN();
	vector<double> d2(d.size(),0);
	vector<int> count(d.size(),0);
	for (int edge=0;edge<=n*2-1;edge++) {
	int beg=0;
	for (int e=0;e<edge;e++) beg+=floor(tree->getDist(e)/precision);
	int end=beg;
	if (edge<n*2-1) end+=floor(tree->getDist(edge)/precision);else end+=floor((timeScale-tree->getNode(n*2-2)->getAge())/precision);
	for (int w=beg;w<end;w++) {
		int a=-5,b=5;
		if (w+a<beg) {b+=beg-w-a;}else
		if (w+b>=end) {a-=w+b-end+1;}
		for (int k=a;k<=b;k++) if (w+k>=beg && w+k<end && w+k<d.size() && w<d2.size()) {d2.at(w)+=d.at(w+k);count.at(w)++;}  //*** There is a danger that something is wrong here because when mapping to a consensus tree sometimes end>d.size....
	}
	}
	for (unsigned int i=0;i<d.size();i++) if (count.at(i)>0) d.at(i)=d2.at(i)/count.at(i);//else cout<<"error in smoothing"<<endl;
}

void DensityOnTree::display(QPainter * painter,vector<double>*x,vector<double>*y,double scale,QBrush brush)
{
  int n=tree->getN();
  for (int i=0;i<n+n-1;i++) {
      int which=0;
      for (int e=0;e<i;e++) which+=floor(tree->getDist(e)/precision);
      for (unsigned int j=0;(i==n*2-2&&j<floor((timeScale-tree->getNode(i)->getAge())/precision))||(i<n*2-2&&j<floor(tree->getDist(i)/precision));j++) {
        double beg=x->at(i)-(j+1)*precision/timeScale;
        painter->fillRect(QRectF(beg,y->at(i)+0.002*((scale<0)?-1:1),precision/timeScale*1.05,d[which]*scale/n),brush);
        which++;
        }
    }
  if (above>0) painter->fillRect(QRectF(0.0,y->at(n+n-2),-0.1,above*scale/n*precision/timeScale*10.0),brush);
}

void DensityOnTree::add2(RecTree *treefrom,int i,bool goingto,double relcon)
{
	vector<int> alle;
	double distto=0;
	if(goingto) alle =treefrom->getAllSampledSeqs(treefrom->getEdge(i)->getEdgeTo());
	else alle =treefrom->getAllSampledSeqs(treefrom->getEdge(i)->getEdgeFrom());
	vector<int> allaffected=tree->getMinEdgeList(alle);

	for(unsigned int j=0;j<allaffected.size();j++) {
	  if(goingto) {
		if(allaffected[j]==tree->getRoot()->getId()) distto=treefrom->getRecEdge(i)->getTimeTo();
		else if(treefrom->getRecEdge(i)->getEdgeTo()==treefrom->getRoot()->getId()) {distto = tree->getNode(allaffected[j])->getDist()*0.99;}
		else distto=treefrom->getRecEdge(i)->getTimeTo()/ treefrom->getNode(treefrom->getRecEdge(i)->getEdgeTo())->getDist() * tree->getNode(allaffected[j])->getDist();
	  }else {
		if(allaffected[j]==tree->getRoot()->getId()) distto=treefrom->getRecEdge(i)->getTimeFrom();
		else if(treefrom->getRecEdge(i)->getEdgeFrom()==treefrom->getRoot()->getId()) {distto = tree->getNode(allaffected[j])->getDist()*0.99;}
		else distto=treefrom->getRecEdge(i)->getTimeFrom()/ treefrom->getNode(treefrom->getRecEdge(i)->getEdgeFrom())->getDist() * tree->getNode(allaffected[j])->getDist();
	  }
	add(allaffected[j],distto,relcon/allaffected.size());
	}
}



//
