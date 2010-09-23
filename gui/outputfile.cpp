#include "outputfile.h"
//
OutputFile::OutputFile(QStringList qstrs,bool makeVectors) 
{
xml.clear();
for (int i=0;i<qstrs.size();i++) {
	if (!qstrs[i].contains("*")&&!qstrs[i].contains("?"))
	file.push_back(new QFile(qstrs[i]));
	else
	{
	if (qstrs[i].lastIndexOf("/")==-1) qstrs[i]="./"+qstrs[i];
	int ind=qstrs[i].lastIndexOf("/");
	QString l=qstrs[i].left(ind);
	QString r=qstrs[i].right(qstrs[i].size()-ind-1);
	QDir dir(l);
	QStringList list=dir.entryList(QStringList(r));
	for (int j=0;j<list.size();j++) {
		cout<<l.toStdString()<<"/"<<list[j].toStdString()<<endl;
		file.push_back(new QFile(l+"/"+list[j]));
	}
	try{
		string ulims=getStdoutFromCommand("ulimit -n");
		int ulim=atoi(ulims.c_str());
		if(ulim<(int)file.size()+1 && ulim>0) {
			cerr<<"Error: You are trying to open more files than your shell allows.  Try running: \"ulimit -n "<<file.size()+1<<"\" from the command line before executing the gui."<<endl;
		}
	}catch(...){cerr<<"Warning: file limit problems.  You may not be able to open all files"<<endl;}
	};
}
outputinitialised=false;
startOver();
//if (makeVectors) countVectors();
}

void OutputFile::reset()
{
outputinitialised=false;
    int currentIteration;
    thetas.clear();
    rhos.clear();
    deltas.clear();
    likelihoods.clear();
    priors.clear();
    numrecedges.clear();
    genorec.clear();
    relgenorec.clear();
    genobeg.clear();
    genoend.clear();
    tmrcas.clear();
    ttotals.clear();
}

void OutputFile::startOver()
{
  currentIteration=-1;
  xml.clear();
  //Open file
  for (unsigned int i=0;i<file.size();i++) {
	file[i]->close();
  if (!file[i]->open(QIODevice::ReadOnly)) {cerr<<"Unable to open file "<<file[i]->fileName().toStdString()<<endl;exit(1);}
  xml.push_back(new QXmlStreamReader(file[i]));
  //Load blocks
  QString str="";
  xml[i]->readNext();
  while (str.compare("Blocks")!=0) {xml[i]->readNext();if (xml[i]->error()!=0) {cerr<<"Invalid file "<<file[i]->fileName().toStdString()<<endl;exit(1);};str=xml[i]->name().toString();}
  xml[i]->readNext();
  QString tmpBlocks=xml[i]->text().toString();
  while( str.compare("Iteration")!=0){
	while (str.compare("comment")!=0 && str.compare("Iteration")!=0  && str.compare("nameMap")!=0  && str.compare("regions")!=0) {xml[i]->readNext();if (xml[i]->error()!=0) {cerr<<"Invalid file "<<file[i]->fileName().toStdString()<<endl;exit(1);};str=xml[i]->name().toString();}
		if(str.compare("comment")==0) {comment=xml[i]->readElementText();str=QString("");}
		if(str.compare("nameMap")==0) {readNames(xml[i]->readElementText());str=QString("");}
		if(str.compare("regions")==0) {addRegions(xml[i]->readElementText());str=QString("");}
	}
 addBlocks(tmpBlocks,i);
 }
// rewind the document
 xml.clear();
 for (unsigned int i=0;i<file.size();i++) {file[i]->close();
  if (!file[i]->open(QIODevice::ReadOnly)) {cerr<<"Unable to open file "<<file[i]->fileName().toStdString()<<endl;exit(1);}
  xml.push_back(new QXmlStreamReader(file[i]));
 }
}

bool OutputFile::getIt(ParamQt * p)
{
   p->setRho(0);p->setTheta(0);p->setLL(0);
   for (unsigned int i=0;i<file.size();i++) {
   int deb=0;
   if (i>0) deb=blocks.split("\n").at(i).split(",").last().toInt();
   xml[i]->readNext();
   while (xml[i]->isStartElement()==false || xml[i]->name().toString().compare("Iteration")!=0) {if (xml[i]->atEnd()) return false;xml[i]->readNext();}
   xml[i]->readNext();
   int start=0,end=0,efrom=0,eto=0;double ato=0,afrom=0;
   while (xml[i]->name().toString().compare("Iteration")!=0) {
  if(xml[i]->error()!=QXmlStreamReader::NoError) {
	cerr<<"XML error of type "<<xml[i]->error()<<":"<<xml[i]->errorString().toStdString()<<endl;
	return true;
  }
    if (xml[i]->isStartElement() && xml[i]->name().toString().compare("Tree")==0) {
      xml[i]->readNext();
      string s=xml[i]->text().toString().toStdString();
      while (s.at(0)==10 || s.at(0)==13) s=s.substr(1,s.length()-1);
      while (s.at(s.size()-1)==10 || s.at(s.size()-1)==13) s=s.substr(0,s.length()-1);
      if (i==0) p->setTreeData(new RecTree(getL(),s,false,false),blocks);
    }
    if (xml[i]->isEndElement()   && xml[i]->name().toString().compare("recedge")==0) {	
	try{p->getTree()->addRecEdge(afrom,ato,start,end,efrom,eto);
	}catch( char * str ){
	cerr<<"Exception from getIt: "<<str<<endl;
	cerr<<"Ignoring this edge!  In Iteration "<<currentIteration<<". Edge details: afrom="<<afrom<<" ato="<<ato<<" efrom="<<efrom<<" eto="<<eto<<" start="<<start<<" end="<<end<<endl;}
    }
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("number")==0) p->setNumber(xml[i]->readElementText().toLong());
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("theta")==0) p->setTheta(p->getTheta()+xml[i]->readElementText().toDouble());
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("delta")==0) p->setDelta(xml[i]->readElementText().toDouble());
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("rho"  )==0) p->setRho  (p->getRho()+xml[i]->readElementText().toDouble());
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("ll"   )==0) p->setLL(p->getLL()+xml[i]->readElementText().toDouble());
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("start")==0) start=deb+xml[i]->readElementText().toInt();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("end"  )==0) end  =deb+xml[i]->readElementText().toInt();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("efrom")==0) efrom=xml[i]->readElementText().toInt();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("eto"  )==0) eto  =xml[i]->readElementText().toInt();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("afrom")==0) afrom=xml[i]->readElementText().toDouble();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().compare("ato"  )==0) ato  =xml[i]->readElementText().toDouble();
    else if (xml[i]->isStartElement() && xml[i]->name().toString().startsWith("est")) {
	p->addConv(xml[i]->name().toString().toStdString(),xml[i]->readElementText().toDouble());
	}
    else if (xml[i]->isStartElement() && xml[i]->name().toString().startsWith("acc")) {
	p->addConv(xml[i]->name().toString().toStdString(),xml[i]->readElementText().toDouble());
	}
    else if (xml[i]->isStartElement() && xml[i]->name().toString().startsWith("prop")) {
	p->addConv(xml[i]->name().toString().toStdString(),xml[i]->readElementText().toDouble());
	}
    xml[i]->readNext();
}// end while loop
   }
   currentIteration++;
   return true;
}

void OutputFile::countVectors()
{
double curNumRec=0.0;
int start=0;
int end=0;
int itercounts=0;
outputinitialised=true;
//
  currentIteration=-1;
  xml.clear();
//Open file
  file[0]->close();
  file[0]->open(QIODevice::ReadOnly);
  xml.push_back(new QXmlStreamReader(file[0]));
  QString str="";
  xml[0]->readNext();

while (!xml[0]->atEnd()) {
  xml[0]->readNext();
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("Blocks" )==0) {int L=xml[0]->readElementText().split(",").last().toInt();genorec=vector<double>(L,0);genobeg=vector<double>(L,0);genoend=vector<double>(L,0);continue;}
  if (xml[0]->isEndElement()&&xml[0]->name().toString().compare("Iteration")==0) {numrecedges.push_back(curNumRec);curNumRec=0.0;continue;}
  if (xml[0]->isEndElement()&&xml[0]->name().toString().compare("recedge"  )==0) {curNumRec++;genobeg[start]++;genoend[end-1]++;for (int i=start;i<end;i++) genorec[i]++;continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("start")==0) {start=xml[0]->readElementText().toInt();continue;} 
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("end"  )==0) {end  =xml[0]->readElementText().toInt();continue;} 
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("theta")==0) {thetas.push_back(xml[0]->readElementText().toDouble());continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("delta")==0) {deltas.push_back(xml[0]->readElementText().toDouble());continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("rho"  )==0) {rhos  .push_back(xml[0]->readElementText().toDouble());continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("ll"   )==0) {likelihoods.push_back(xml[0]->readElementText().toDouble());continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("prior")==0) {priors.push_back(xml[0]->readElementText().toDouble());continue;}
  if (xml[0]->isStartElement()&&xml[0]->name().toString().compare("Tree" )==0) {itercounts++; string s=xml[0]->readElementText().toStdString();while (s.at(0)==10 || s.at(0)==13) s=s.substr(1,s.length()-1);while (s.at(s.size()-1)==10 || s.at(s.size()-1)==13) s=s.substr(0,s.length()-1);Tree * t=new Tree(s,false);tmrcas.push_back(t->getNode(t->getN()*2-2)->getAge());ttotals.push_back(t->getTTotal());delete(t);continue;}
}
	for(unsigned int i=0;i<genorec.size();i++) {
		genorec[i]/=itercounts;
		genobeg[i]/=itercounts;
		genoend[i]/=itercounts;
	}
startOver();
}

vector<double>* OutputFile::getRhoOverTheta()
{
vector<double>*v=new vector<double>();
for (unsigned int i=0;i<rhos.size();i++) v->push_back(rhos[i]/thetas[i]);
return v;
}

vector<double>* OutputFile::getRhoPerSite()
{
int L=getL();
int b=getB();
vector<double>*v=new vector<double>();
for (unsigned int i=0;i<rhos.size();i++) v->push_back(rhos[i]/(deltas[i]*b+L-b));
return v;
}

vector<double>* OutputFile::getThetaPerSite()
{
int L=getL();
vector<double>*v=new vector<double>();
for (unsigned int i=0;i<rhos.size();i++) v->push_back(thetas[i]/L);
return v;
}

vector<double>* OutputFile::getPosteriors()
{
vector<double>*v=new vector<double>();
for (unsigned int i=0;i<likelihoods.size();i++) v->push_back(likelihoods[i]+priors[i]);
return v;
}

vector<double>* OutputFile::getRoverM(ParamQt*param)
{
vector<double>*v=new vector<double>();
startOver();
while (getIt(param)) {
v->push_back(param->getRM());
}
startOver();
return v;

/*int L=getL();
int b=getB();
vector<double>*v=new vector<double>();
//for (unsigned int i=0;i<rhos.size();i++) v->push_back(rhos[i]/thetas[i]*L*deltas[i]/(deltas[i]*b+L-b)*0.75*(1.0-exp(-4.0*thetas[i]/L)));
for (unsigned int i=0;i<rhos.size();i++) v->push_back(rhos[i]/thetas[i]*L*deltas[i]/(deltas[i]*b+L-b)*3.0*thetas[i]/(3.0*L+4.0*thetas[i]));
return v;*/
}

void OutputFile::makeCF(ParamQt*param)
{
	vector<vector<double> > *v=new vector<vector<double> >(0,vector<double>(0.0));
	startOver();
	int count=0;
	while (getIt(param)) {
		param->makeCF(v);
		count++;
	}
	param->setCF(v,count);
	startOver();
}

vector<double>* OutputFile::getGenoRec(int id,bool getto) {
vector<double> * res=new vector<double>(genorec.size(),0);
QFile f(file[0]->fileName());
f.open(QIODevice::ReadOnly);
QXmlStreamReader x(&f);
int start=0;
int end=0;
int edge=0;int efrom=0;
while (!x.atEnd()) {
  x.readNext();
  if (x.isEndElement()&&x.name().toString().compare("recedge"  )==0) {
	if (edge==id && getto) for (int i=start;i<end;i++) (res->at(i))++;
	else if(efrom==id && !getto) for (int i=start;i<end;i++) (res->at(i))++;
  	continue;
  }
  if (x.isStartElement()&&x.name().toString().compare("start")==0) {start=x.readElementText().toInt();continue;} 
  if (x.isStartElement()&&x.name().toString().compare("end"  )==0) {end=x.readElementText().toInt();continue;} 
  if (x.isStartElement()&&x.name().toString().compare("efrom")==0) {efrom=x.readElementText().toInt();continue;} 
  if (x.isStartElement()&&x.name().toString().compare("eto")==0) {edge=x.readElementText().toInt();continue;} 
}
f.close();
for (unsigned int i=0;i<res->size();i++) {res->at(i)/=thetas.size();}
return res;
}

vector<double>* OutputFile::getRelGenoRec(ParamQt*param,int id) {
vector<double> * res=new vector<double>(genorec.size(),0);
QFile f(file[0]->fileName());
f.open(QIODevice::ReadOnly);
QXmlStreamReader x(&f);
double treedist;
startOver();
int counts=0;
while (getIt(param)) {
	for(int i=0;i<(int)param->getRecTree()->numRecEdge();i++){
		if(param->getRecTree()->getEdge(i)->getEdgeTo()==id){
		treedist=param->getRecTree()->getEdgeTreeTime(i);
		for (unsigned int j=param->getRecTree()->getEdge(i)->getStart();j<param->getRecTree()->getEdge(i)->getEnd();j++) (res->at(j))+=treedist/param->getRecTree()->getTTotal();
		counts++;
		}
	}
}
f.close();
for (unsigned int i=0;i<res->size();i++) res->at(i)/=counts;
return res;
}

vector<double>* OutputFile::getRelGenoRec(ParamQt*param) {
vector<double> * res=new vector<double>(genorec.size(),0);
QFile f(file[0]->fileName());
f.open(QIODevice::ReadOnly);
QXmlStreamReader x(&f);
double treedist;
startOver();
int counts=0;
while (getIt(param)) {
	for(int i=0;i<(int)param->getRecTree()->numRecEdge();i++){
		treedist=param->getRecTree()->getEdgeTreeTime(i);
		for (unsigned int j=param->getRecTree()->getEdge(i)->getStart();j<param->getRecTree()->getEdge(i)->getEnd();j++) (res->at(j))+=treedist/param->getRecTree()->getTTotal();
		counts++;
	}
}
f.close();
for (unsigned int i=0;i<res->size();i++) res->at(i)*=param->getDelta()/counts;
return res;
}


void OutputFile::readNames(QString str){
	names.clear();
	QStringList list1 = str.split(";");
	for(unsigned int i=0;i<list1.size();i++){
		QStringList list2 = list1[i].split(",");
		if(list2.size()>1){
			int index=list2[0].toInt();
			names<<list2[1];
		}
	}
}


string OutputFile::getStdoutFromCommand(string cmd)
{
  // setup
  string data;
  FILE *stream;
  int MAX_BUFFER=10000;
  char buffer[MAX_BUFFER];

  // do it
  stream = popen(cmd.c_str(), "r");
  while ( fgets(buffer, MAX_BUFFER, stream) != NULL )
    data.append(buffer);
  pclose(stream);

  // exit
  return (data);
}
