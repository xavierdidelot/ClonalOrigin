#ifndef GELMANRUBINIMPL_H
#define GELMANRUBINIMPL_H
//
#include "ui_gelmanrubin.h"
#include <QtXml>
#include <iostream>
#include <vector>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include "paramqt.h"
#include "../warg/src/tree.h"
#include "../warg/src/rectree.h"
#include "../warg/src/param.h"
#include "../warg/src/data.h"
#include "outputfile.h"
#include "paramtreecons.h"
//

using namespace std;

class GelmanRubinImpl : public QDialog, public Ui::GelmanRubin
{
Q_OBJECT
public:
	GelmanRubinImpl( QWidget * parent = 0, Qt::WFlags f = 0 );
	void compute(ParamQt*param,OutputFile*outputFile,QStringList*others, ostream* out=NULL,bool getparams=true,bool getnumedges=true,bool getpairwisedists=false);
	void computeTree(ParamQt*,QStringList*,ostream* out=NULL);
	void outputTracer(ParamQt*param,OutputFile*outputfile,QString*qstr,bool csv=false,bool getparams=true,bool getnumedges=true,bool getpairwisedists=false);
	static QString compareTrueTree(ParamQt*,QString,QString);

	inline void setFiles(ParamQt*param,OutputFile*outputFile,QStringList others){
		this->param=param;
		this->outputFile=outputFile;
		this->others=others;
	}
	inline void showOptions(){
		groupBox->setHidden(false);
	}
	inline void setExport(){
		groupBox->setHidden(false);
		table->setHidden(true);
		exportButton->setHidden(true);
		label->setHidden(true);
		pushGo->setHidden(true);
		doGR=false;
		this->setFixedHeight(170);
		this->setWindowTitle(QString("Export"));
	}
	inline void setSep(string s){sep=s;}
protected:
	void nameTable(ParamQt*p,bool getparams,bool getnumedges,bool getpairwisedists);
	vector<double> extractInfo(ParamQt*p, bool csv=false,bool getparams=true,bool getnumedges=true,bool getpairwisedists=false);
	double test(vector< vector<double> >*data);
	ParamQt*param;
	OutputFile*outputFile;
	QStringList others;
	ostream* out;
	string sep;
	bool doGR;
private slots:
	void on_exportButton_clicked();
	void on_buttonBox_accepted();
	void on_buttonBox_rejected();
	void on_pushGo_clicked();

};


#endif
