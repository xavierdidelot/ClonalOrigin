#ifndef PARAMMR_H
#define PARAMMR_H
//
#include "paramcons.h"
//
class ParamMR : public ParamQt
{

public:
	ParamMR();
	virtual ~ParamMR();
	void account();
	void consensus(int cutoff,int site=0);
	void consensusAllSites(int cutoff);
	void makeDensity(ParamCons*paramcons);
	void toCSV(QTextStream*out,int cutoff,int step);
	void toArtemis(QTextStream*out,int cutoff);
	QString toString();
    void correctForPrior();
protected:
	int its;
	int l;
	int coef;
	vector<vector<RecEdge> > edges;
	vector<RecEdge> means;
	bool sameEdge(RecEdge*,RecEdge*);
	double dist(int,double,int,double);
	RecEdge average(RecEdge*r1,RecEdge*r2,int coef);
	void gotowards(int e1,double a1,int e2,double a2,int*e3,double*a3,double dist);
    vector<vector<vector<int> > > states;
};
#endif
