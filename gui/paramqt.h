#ifndef PARAMQT_H
#define PARAMQT_H
//
#include <QtXml>
#include <QtGui>
#include "gsl/gsl_math.h"
#include "../warg/src/param.h"
//
using namespace weakarg;
using namespace std;

class ParamQt :public QObject, public Param {
      Q_OBJECT

public:
	ParamQt();
	~ParamQt();
	virtual void display(QPaintDevice * qpd); 
	void getXY(vector<double>**xx,vector<double>**yy);
	void makeYforLeaves(vector<double> * v,Node * cur,int *x);
	void displayTree(QPainter *painter,vector<double>*x,vector<double>*y);
	void displayEdges(QPainter *painter,vector<double>*x,vector<double>*y);
	inline void setGene(int g) {gene=g;}
	inline int getGene() {return gene;}
	void setBlocks(QString qstr);
    inline void setTreeData(RecTree * t,QString strblocks) {
        delete(rectree);
        rectree=t;
        if (data==NULL) {
	vector<int> blocks;
	QStringList qsl=strblocks.split(",");
	for (int i=0;i<qsl.size();i++) blocks.push_back(qsl.at(i).toInt());
	data=new Data(rectree->getN(),blocks);}
    }///<Sets the tree
	inline void clearTreeData(){
		if(data==NULL) delete(data);
	}
	bool isCons;
	virtual inline void incrTimeScale() {timeScale*=1.1;};
	virtual inline void decrTimeScale() {timeScale/=1.1;};
	inline void setTimeScale(double ts) {timeScale=ts;};
	inline double getRateScale() {return rateScale;}
	inline void setRateScale(double r) {rateScale=r;} 
	inline double getTimeScale() {return timeScale;}
	double getRM();
	//inline void clearConv() {convnames.clear();convdata.clear();}
	inline void addConv(string name, double data){
		for(unsigned int i=0;i<convnames.size();i++) {
			if(name.compare(convnames.at(i))==0) {convdata[i]=data;return;}
		}
		convnames.push_back(name);
		convdata.push_back(data);
	}
	inline string getConvName(int index){
		if(index>=(int)convnames.size()){cerr<<"Error in paramqt: index "<<index<<" doesn't exist in convnames"<<endl;throw;};
		return convnames.at(index);}
	inline double getConvData(int index){
		if(index>=(int)convdata.size()){cerr<<"Error in paramqt: index "<<index<<" doesn't exist in convdata"<<endl;throw;};
		return convdata.at(index);}
	inline int countConv(){return convdata.size();}
	int getIdAt(int x,int y,QPaintDevice * qpd);
	inline void setDisplayTree(RecTree *intree,bool isnew=false){
		displaytree=intree;
		if(isnew) displayset=true;
	}///*Sets the tree that will be displayed
	inline void newDisplayTree(RecTree *intree,bool forceages=false){
		if(displayset) delete(displaytree);
		//displaytree=new RecTree(intree,intree->newick(64),false,forceages);
		displaytree=new RecTree(*intree);
		displayset=true;
	}///*Sets the tree that will be displayed to a new tree copied from intree
	inline RecTree * getDisplayTree(){
		return(displaytree);
	}///*Gets the tree that will be displayed
	inline void unsetDisplayTree(){
		if(displayset) delete(displaytree);
		displayset=false;
	}///*removes the display tree
	inline bool displaySet(){
		return(displayset);
	}///*returns whether the tree is set
	inline void setNumber(long i){
		iteration=i;
	}///* Sets the iteration we are on
	inline long getNumber(){
		return(iteration);
	}///* Sets the iteration we are on
	inline void toggleRecView(){
		recview++;
		if(recview>=3) recview=0;
	}
	inline void setNameType(int i){
		nametype=i;
	}
	inline void setRecView(int i){if(i>=0 &&i<3) recview=i;}	
	inline int getRecView(){return(recview);}
	RecTree *displaytree;
	int displayset;
	double rateScale;
	double timeScale;
	long iteration;
	int recview;
	int labview;
	int nametype;
	vector<double> convdata;// convergence diagnostics
	vector<string> convnames;// convergence diagnostics
    	QStringList labels;
    	QStringList names;
	inline QString isolateName(int i){
		if(i>=names.size() || i<0 || nametype<0) { return(QString::number(i));
		}else {
		if(nametype==1) {
			QStringList tmp=names[i].split("+");
			tmp.erase(tmp.begin());
			return(tmp.join(" "));
		}else if(nametype==2) {
			QStringList tmp=names[i].split("+");
			tmp.erase(tmp.begin());
			tmp=tmp.join(" ").split(".");
			tmp.erase(tmp.begin()+tmp.size()-1);
			return(tmp.join(" "));
		}

		return(names[i]);
		}
    	}
	inline QString nodeLabel(int i){
		if(i>=labels.size() || i<0) { return(QString(""));
		}else return(labels[i]);
    	}
	void setNames(QStringList qn){
		names=qn;
	}
	void setLabels(QStringList qn){
		labels=qn;
	}
	void makeCF(vector <vector<double> > *v);///<Accounts for clonal frame proportions in v
	void setCF(vector <vector<double> > *v,int count);///< Sets the ClonalFrame proportions as labels
	int gene;//Which gene to show; -1 for all
	int lastCommonAncestor(int s1, int s2);///< Returns the last common ancestor between two individuals
	vector<int>  consistentAgeList(vector<double> *res);///< Returns a list of node ordersthat is consistent if -f option is used, and puts their ages in res
	double pairwiseDistance(int s1, int s2);///<Returns the pairwise distance between sequence s1 and s2
	vector<double> pairwiseDistanceList();///<List of the unique pairwise distances i9n order for(c1=0..N-1 for(c2=c1+1..N-1)).
	vector<vector<double> > pairwiseDistanceMatrix(bool print=false);///< returns a pairwise distance matrix for the current tree
	int recCount(int s1, int s2);///< Gets the number of recombination events from s2 to s1 (or ancestors) before they coalesce
	vector<vector<int> > recCountMatrix(bool print);///< Gets the pairwise count of recombination events
	vector<vector<double> > recPriorMatrix();///< recombination prior matrix, pairwise 
};
#endif
