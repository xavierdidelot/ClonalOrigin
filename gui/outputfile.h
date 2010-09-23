#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H
#include <QtXml>
#include <QtGui>
#include "paramqt.h"
//
using namespace std;
//
class OutputFile  
{
public:
	OutputFile(QStringList qstrs,bool makeVectors=true);
	inline int getL() {return blocks.split(",").last().toInt();}
	inline int getB() {return blocks.split(",").size()-1;}
	inline QString getBlocks() {return blocks;}
	bool getIt(ParamQt * p);
	void startOver();
	void reset();
	inline int getCurIt() {return currentIteration;}
    void addOtherData(string name,double val);
    inline vector<double>*getThetas(){vector<double> *v=new vector<double>(thetas);return v;}
    inline vector<double>*getRhos(){vector<double> *v=new vector<double>(rhos);return v;}
    inline vector<double>*getDeltas(){vector<double> *v=new vector<double>(deltas);return v;}
    inline vector<double>*getLikelihoods(){vector<double> *v=new vector<double>(likelihoods);return v;}
    inline vector<double>*getPriors(){vector<double> *v=new vector<double>(priors);return v;}
    inline vector<double>*getNumRecEdges(){vector<double> *v=new vector<double>(numrecedges);return v;}
    inline vector<double>*getGenoRec(){vector<double> *v=new vector<double>(genorec);return v;}
    vector<double>*getRelGenoRec(ParamQt*param);
    vector<double>*getGenoRec(int id,bool getto);
    vector<double>*getRelGenoRec(ParamQt*param,int id);
    inline vector<double>*getGenoBeg(){vector<double> *v=new vector<double>(genobeg);return v;}
    inline vector<double>*getGenoEnd(){vector<double> *v=new vector<double>(genoend);return v;}
    vector<double>*getRhoOverTheta();
    vector<double>*getRoverM(ParamQt*param);
    vector<double>*getRhoPerSite();
    vector<double>*getThetaPerSite();
    vector<double>*getPosteriors();
    vector<double>*getTMRCA() {vector<double> *v=new vector<double>(tmrcas);return v;}
    vector<double>*getTTotal() {vector<double> *v=new vector<double>(ttotals);return v;}
    inline bool isinitialised(){return(outputinitialised);}
    void countVectors();
    inline QString getFileName() {return file[0]->fileName();}
    inline QString getComment(){return(comment);}
    void readNames(QString str);
    inline QStringList getNames(){return(names);};
    void makeCF(ParamQt*param);
    inline void addRegions(QString str){
	QStringList strl=str.split(",");
	for(unsigned int i=0;i<strl.size();i++) regions.push_back(strl[i].toInt());
    }// adds to the list of regions, keeping the list sorted and unique
    inline vector<int> getRegions(){
	return(regions);
    }
   inline void addBlocks(QString tmpBlocks,int i){
	previousL.push_back(getL());
	if (i==0) blocks=tmpBlocks; else {
	int L=getL();
	QStringList qstrl=tmpBlocks.split(",");
	for (int j=1;j<qstrl.size();j++) blocks.append(",").append(QString::number(L+qstrl[j].toInt()));
	blocks.append("\n");
	};
    }

protected:
    bool outputinitialised;
    QString comment;
    vector<QFile*>file;
    vector<QXmlStreamReader*> xml;
    QString blocks;
    int currentIteration;
    vector<int> previousL;
    vector<double>thetas;
    vector<double>rhos;
    vector<double>deltas;
    vector<double>likelihoods;
    vector<double>priors;
    vector<double>numrecedges;
    vector<double>genorec;
    vector<double>relgenorec;
    vector<double>genobeg;
    vector<double>genoend;
    vector<double>tmrcas;
    vector<double>ttotals;
    vector<int> regions;
    QStringList names;
    std::string getStdoutFromCommand(std::string cmd);
};
#endif
