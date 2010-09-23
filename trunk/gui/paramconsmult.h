#ifndef PARAMCONSMULT_H
#define PARAMCONSMULT_H
//
#include "paramqt.h"
#include "densityontree.h"
//
class ParamConsMult  : public ParamQt
{
public:
	ParamConsMult(QStringList*nodes,QStringList*colors,bool denDep,bool colDep);
	virtual ~ParamConsMult();
	void account();
	void display(QPaintDevice * qpd);
	inline bool accept(int node1,int node2) {if (node1==node2) return true;
	if (node1==2*rectree->getN()-2) return false;
	return accept(rectree->getNode(node1)->getFather()->getId(),node2);};
	inline void set() {isSet=true;for (unsigned int i=0;i<dps.size();i++) {dps[i]->smooth();dms[i]->smooth();}}
	inline void incrTimeScale() {QMessageBox::about(0, "Information","Can't change scale while looking at rates.");};
	inline void decrTimeScale() {QMessageBox::about(0, "Information","Can't change scale while looking at rates.");}; 
protected:
	vector<int> nodes;
	vector<QBrush> brushes;
	int its;
	bool isSet;
	bool denDep;
	bool colDep;
	vector<DensityOnTree*> dps;///<Density of gain
	vector<DensityOnTree*> dms;///<Density of loss
};
#endif
