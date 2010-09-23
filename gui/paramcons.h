#ifndef PARAMCONS_H
#define PARAMCONS_H
//
#include "paramqt.h"
#include "densityontree.h"
//
class ParamCons  : public ParamQt
{
public:
	ParamCons();
	virtual ~ParamCons();
	void account(int id=-1,bool getto=true);
	void reset();
	void accountrelative(int id=-1,bool getto=true);
	void display(QPaintDevice * qpd);

	inline void set() 
	{
	isSet=true;
	dp->smooth();dm->smooth();
	}
	inline void incrTimeScale() {QMessageBox::about(0, "Information","Can't change scale while looking at rates.");};
	inline void decrTimeScale() {QMessageBox::about(0, "Information","Can't change scale while looking at rates.");}; 
protected:
	int its;
	bool isSet;
	DensityOnTree * dp;///<Density of gain
	DensityOnTree * dm;///<Density of loss
};
#endif
