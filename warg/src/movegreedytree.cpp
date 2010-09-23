#include "movegreedytree.h"

//#define DEBUG

using namespace std;
namespace weakarg
{

MoveGreedyTree::MoveGreedyTree(Param * p,double a)
        : Move(p,a)
{
    description="Greedily update the tree";
    name="GreedyTree";
}


int MoveGreedyTree::recCount(int s1, int s2)
{
	RecTree * rectree=param->getRecTree();
	vector<int> r1,r2;
	int lcaindexr1=-1,lcaindexr2=-1;/// index of last common ancestor in list from s1
	int reccounted=0;
// get the list of edges to the MRCA
	r1.push_back(s1);
	while(r1.back()!=rectree->getRoot()->getId()) {
		r1.push_back(rectree->getNode(r1.back())->getFather()->getId());
	}
	r2.push_back(s2);
	while(r2.back()!=rectree->getRoot()->getId()) {
		r2.push_back(rectree->getNode(r2.back())->getFather()->getId());
	}
	for(unsigned int i=1;i<=max(r1.size(),r2.size());i++) {
		if(r1[r1.size()-i]==r2[r2.size()-i]) {
			lcaindexr1=r1.size()-i;
			lcaindexr2=r2.size()-i;
		}else break;
	}
// get the recombination affecting those edges
	for(long c1=0;c1<rectree->numRecEdge();c1++) {
		for(int c2=0;c2<lcaindexr1;c2++) {
		  if(rectree->getEdge(c1)->getEdgeTo()==r1[c2] || rectree->getEdge(c1)->getEdgeTo()==s1) {
			for(int c3=0;c3<lcaindexr2;c3++) {
			  if(rectree->getEdge(c1)->getEdgeFrom()==r2[c3] || rectree->getEdge(c1)->getEdgeFrom()==s2) {reccounted++;c3=lcaindexr2;c2=lcaindexr1;}
			}
		  }
		}
	}
	return(reccounted);
}


vector<int> MoveGreedyTree::edgesBetween(int s1, int s2)
{
	RecTree * t=param->getRecTree();
	vector<int> r1,r2;
	int lcaindexr1=-1,lcaindexr2=-1;/// index of last common ancestor in list from s1
// get the list of edges to the MRCA
	r1.push_back(s1);
	while(r1.back()!=t->getRoot()->getId()) {
		r1.push_back(t->getNode(r1.back())->getFather()->getId());
	}
	r2.push_back(s2);
	while(r2.back()!=t->getRoot()->getId()) {
		r2.push_back(t->getNode(r2.back())->getFather()->getId());
	}
	for(unsigned int i=1;i<=max(r1.size(),r2.size());i++) {
		if(r1[r1.size()-i]==r2[r2.size()-i]) {
			lcaindexr1=r1.size()-i;
			lcaindexr2=r2.size()-i;
		}else break;
	}
	vector<int> ret;
	for(int c1=0;c1<lcaindexr1;c1++) {
		ret.push_back(r1[c1]);
	}
	for(int c1=0;c1<lcaindexr2;c1++) {
		ret.push_back(r2[c1]);
	}
	return(ret);
}

int MoveGreedyTree::move()
{
	vector<double> dists=calcDists();
// Now we know the pairwise mutation counts, we just need to update the tree.  First, calculate the new ages
	vector<double> newage=calcAges(dists);
// Apply the changes to the tree
	applyChanges(newage);

	return(1);
}

vector<vector<bool> >* MoveGreedyTree::localClonalFrame()
{
	RecTree * t=param->getRecTree();
	long L = t->getL();
	int N = t->getN();
	int npairs = N * (N -1)/2;
	vector<vector<bool> > * lcf = new vector <vector<bool> > (npairs,vector<bool>(L,true)) ;
	vector<vector <int> > slist;
// create a list of which edges are between each pair of sequences
	for(int c1=0;c1<N;c1++) for(int c2=c1+1;c2<N;c2++) slist.push_back(edgesBetween(c1,c2));
// go through every edge	
	for(long c1=0;c1<t->numRecEdge();c1++){
	  if(t->getRecEdge(c1)->getEdgeTo() != t->getRecEdge(c1)->getEdgeFrom()) {
// go through every pair of sequences
	  for(unsigned int c2=0;c2<slist.size();c2++){for(unsigned int c3=0;c3<slist[c2].size();c3++){
// check if the edge affects this pair of sequences
	    if(t->getRecEdge(c1)->getEdgeTo()==slist[c2][c3]) {
// remove the affected sites from the local clonal frame
		for(unsigned long c4=t->getRecEdge(c1)->getStart();c4<t->getRecEdge(c1)->getEnd();c4++) lcf->at(c2)[c4]=false;
	    }
	  }}
	  }
	}
	return(lcf);
}

 vector<double> MoveGreedyTree::calcDists(vector<double> nmuts,vector<double> nsites)
{
	long L = param->getRecTree()->getL();
	vector<double> dists(nmuts.size(),0.0);
	for(unsigned int c1=0;c1<dists.size();c1++) { 
		nsites[c1]+=0.1;nmuts[c1]+=0.001;
		dists[c1]=((double)nmuts[c1])/nsites[c1]*L/param->getTheta() 
 +gsl_ran_gaussian(rng,0.01 * (double)nmuts[c1]/nsites[c1]*L/param->getTheta());
		if(dists[c1]<0) dists[c1]= -dists[c1];
// Its necessary to perturb the distances by a small amount to prevent clashes
	}
	return(dists);
}

vector<double> MoveGreedyTree::calcDists(double * esttheta)
{
	RecTree * t=param->getRecTree();
	int N = t->getN();

	vector<vector<bool> >* lcf = localClonalFrame();
	vector<double> nsites(lcf->size(),0.0);
	vector<double> nmuts(lcf->size(),0.0);
	vector<int> nrecs(lcf->size(),0);
	long pairon=0;
	for(int c1=0;c1<N;c1++){ for(int c2=c1+1;c2<N;c2++){
	  nrecs[pairon]=recCount(c1,c2);
	  for(unsigned int c3=0;c3<lcf->at(pairon).size();c3++) {
	    if(lcf->at(pairon)[c3]){
		nsites[pairon]+=1.0;
		if(param->getData()->get(c1,c3)!=param->getData()->get(c2,c3))nmuts[pairon]+=1.0;
	    }
	  }
	  pairon++;
	}}

	vector<double> dists = calcDists(nmuts,nsites);
	if(esttheta!=NULL){
		//*esttheta = empiricalTheta(nsites,nmuts,dists);
	}
	delete(lcf);
	return(dists);
}
/*
double MoveGreedyTree::empiricalTheta(vector<double> nsites,vector<double> nmuts,vector<double> dists)
{
	RecTree * t=param->getRecTree();
	int N = t->getN();

	vector<double> newage(2*N - 1,0.0);
	for(int c1=N;c1<2*N - 1;c1++){
	   vector<int> leftchildren=getAllTips(t->getNode(c1)->getLeft()->getId(),N,t);
	   vector<int> rightchildren=getAllTips(t->getNode(c1)->getRight()->getId(),N,t);

// calculate the average mutation distance between the two children
	  double sumdist=0;
	  int paircount=0;
	  for(unsigned int c2=0;c2<leftchildren.size();c2++) {
	    for(unsigned int c3=0;c3<rightchildren.size();c3++) {
		int pind=getPairIndex(leftchildren[c2],rightchildren[c3],N);
		  sumdist+=dists[pind];
		  paircount++;
	    }
	  }
	  if(sumdist>0) {newage[c1]=sumdist/paircount;
	  }else newage[c1]=t->getNode(c1)->getAge();
	}
	return(newage);

}
*/
vector<double> MoveGreedyTree::calcAges(vector<double> dists)
{
	RecTree * t=param->getRecTree();
	int N = t->getN();

	vector<double> newage(2*N - 1,0.0);
	for(int c1=N;c1<2*N - 1;c1++){
	   vector<int> leftchildren=getAllTips(t->getNode(c1)->getLeft()->getId(),N,t);
	   vector<int> rightchildren=getAllTips(t->getNode(c1)->getRight()->getId(),N,t);

// calculate the average mutation distance between the two children
	  double sumdist=0;
	  int paircount=0;
	  for(unsigned int c2=0;c2<leftchildren.size();c2++) {
	    for(unsigned int c3=0;c3<rightchildren.size();c3++) {
		int pind=getPairIndex(leftchildren[c2],rightchildren[c3],N);
		  sumdist+=dists[pind];
		  paircount++;
	    }
	  }
	  if(sumdist>0) {newage[c1]=sumdist/paircount;
	  }else newage[c1]=t->getNode(c1)->getAge();
	}
	return(newage);
}


void MoveGreedyTree::applyChanges(vector<double> newage)
{
	RecTree * t=param->getRecTree();
	int N = t->getN();
	vector<int> map(2*N - 1,0);
	if(map.size()!=newage.size()) {cerr<<"MoveGreedyTree::applyChanges Error: List of ages of size "<<newage.size()<<" but need size "<<map.size()<<endl; throw("MoveGreedyTree Error: Age size incorrect");}
	for(int c1=N;c1<2*N - 1;c1++) map[c1]=c1;
//cout<<"Newage : oldage"<<endl; 
//for(int c1=0;c1<newage.size();c1++)cout<<newage[c1]<<" : "<<t->getNode(c1)->getAge()<<endl;

	for(int c1=N;c1<2*N - 1;c1++){
	  if (newage[c1]<0.000000000001)newage[c1]=0.000000000001;
	  vector<int> samplespace;
    	  vector<int> listLR;
	  int whichto=-1,tmp;
	  double dist = newage[c1] - t->getNode(map[c1])->getAge();
	  if(dist>0.001 || dist < -0.001) {
	    t->chooseNodePath(map[c1], dist, &listLR, &samplespace);
	    tmp=t->moveNodeTime(map[c1],&whichto,dist,-1,&listLR,&samplespace);
	    if(map[c1]!=whichto) {
		int olde=map[c1];
		map.erase(map.begin() + c1);
		map.insert(map.begin() + whichto,olde);
	    }
	  }
	}
}

MoveGreedyTree::~MoveGreedyTree()
{}



} // end namespace weakarg

