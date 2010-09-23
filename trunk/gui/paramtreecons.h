#ifndef PARAMTREECONS_H
#define PARAMTREECONS_H
//
#include "paramqt.h"
//
class ParamTreeCons : public ParamQt
{
class Cluster{public:int i;double d;Cluster(){};Cluster(int a,double b) {i=a;d=b;}};
public:
	ParamTreeCons();
	virtual ~ParamTreeCons();
	void account();
	void consensus(int cutoff);
	void consensusExt();
	static void makeKey(Node*,QString*);
protected:
	QHash<QString, Cluster> hash;
	int its;
	QString buildsubtree(QString,vector<QString>*);
	bool compat(QString,QString);
};
#endif
