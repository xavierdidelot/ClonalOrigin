#include "parammr.h"
//
ParamMR::ParamMR() 
{
	its=0;
	coef=1;
}

ParamMR::~ParamMR()
{
}

void ParamMR::account()
{
    if (data->getL()>100000) coef=10;
    if (its==0) {int n=data->getN()*2-1;l=data->getL()/coef;states=vector<vector<vector<int> > >(n,vector<vector<int> >(n,vector<int>(l,0)));}
    its++;
    for (int j=0;j<rectree->numRecEdge();j++) {
    int    start=rectree->getRecEdge(j)->getStart()/coef;
    int    end  =rectree->getRecEdge(j)->getEnd  ()/coef;
    int    efrom=rectree->getRecEdge(j)->getEdgeFrom();
    int    eto  =rectree->getRecEdge(j)->getEdgeTo  ();
    for (int site=start;site<end;site++) states[efrom][eto][site]++;
   }
}

void ParamMR::consensus(int cutoff,int site)
{
while (rectree->numRecEdge()>0) rectree->remRecEdge(0);
int n=data->getN()*2-1;
for (int i=0;i<n;i++) for (int j=0;j<n;j++){
if (states[i][j][site]>=cutoff*its/100)
	try{rectree->addRecEdge(rectree->getNode(i)->getDist()/2+0.1*(i==n-1),rectree->getNode(j)->getDist()/2,site,site,i,j);
	}catch( char * str ){
	cerr<<"Exception from loadMR: "<<str<<endl;
	cerr<<"Ignoring this edge!"<<endl;} 
if (i>=data->getN() && states[i][j][site]<cutoff*its/100 && 
	states[rectree->getNode(i)->getLeft ()->getId()][j][site]<cutoff*its/100 && 
	states[rectree->getNode(i)->getRight()->getId()][j][site]<cutoff*its/100 && 
	states[i][j][site]+states[rectree->getNode(i)->getLeft()->getId()][j][site]+states[rectree->getNode(i)->getRight()->getId()][j][site]>=cutoff*its/100)
	try{rectree->addRecEdge(0,rectree->getNode(j)->getDist()/2,site,site,i,j);
	}catch( char * str ){
	cerr<<"Exception from loadMR: "<<str<<endl;
	cerr<<"Ignoring this edge!"<<endl;}
}
}

void ParamMR::consensusAllSites(int cutoff)
{
while (rectree->numRecEdge()>0) rectree->remRecEdge(0);
int n=data->getN()*2-1;
bool cond1,cond2;//Conditions to initiate and to continue a recedge
int start=-1;//First time we had cond1
int end=-1;//Last time we had cond2
int end2=-1;//Last time we had cond1
int c1=cutoff*its/100;
int c2=(cutoff-10)*its/100;
for (int i=0;i<n;i++) for (int j=0;j<n;j++) for (int mode=0;mode<2;mode++) for (int site=0;site<=l;site++) {
if (start!=-1 && site==l) end=site-1;
if (end!=-1) {
	if (end-start>20)
	{if (mode==0) rectree->addRecEdge(rectree->getNode(i)->getDist()/2+0.1*(i==n-1),rectree->getNode(j)->getDist()/2,start,end,i,j);
	else rectree->addRecEdge(0,rectree->getNode(j)->getDist()/2,start,end,i,j);}
	start=-1;
	end=-1;
}
if (site==l) continue;
if (mode==0) {
cond1=states[i][j][site]>=c1;
cond2=states[i][j][site]>=c2;}
else {
cond1=i>=data->getN() && states[i][j][site]<c2 && 
	states[rectree->getNode(i)->getLeft ()->getId()][j][site]<c2 && 
	states[rectree->getNode(i)->getRight()->getId()][j][site]<c2 && 
	states[i][j][site]+states[rectree->getNode(i)->getLeft()->getId()][j][site]+states[rectree->getNode(i)->getRight()->getId()][j][site]>=c1;
cond2=i>=data->getN() && states[i][j][site]<c2 && 
	states[rectree->getNode(i)->getLeft ()->getId()][j][site]<c2 && 
	states[rectree->getNode(i)->getRight()->getId()][j][site]<c2 && 
	states[i][j][site]+states[rectree->getNode(i)->getLeft()->getId()][j][site]+states[rectree->getNode(i)->getRight()->getId()][j][site]>=c2;
}
if (cond1)
	{
		if (start==-1) start=site;
		end2=site;
		continue;
	}
if (start!=-1 && !cond2) end=end2;
}
}

bool ParamMR::sameEdge(RecEdge*r1,RecEdge*r2)
{
	if (r1->getStart()>=r2->getEnd()) return false;
	if (r2->getStart()>=r1->getEnd()) return false;
	if (dist(r1->getEdgeFrom(),r1->getTimeFrom(),r2->getEdgeFrom(),r2->getTimeFrom())>rectree->getTTotal()/50.0) return false;
	if (dist(r1->getEdgeTo  (),r1->getTimeTo  (),r2->getEdgeTo  (),r2->getTimeTo  ())>rectree->getTTotal()/50.0) return false;
	return true;
}

double ParamMR::dist(int e1,double a1,int e2,double a2){
	if (e1==e2) return fabs(a1-a2);
	if (rectree->getNode(e1)->getAge()+a1>rectree->getNode(e2)->getAge()+a2) {swap(e1,e2);swap(a1,a2);}
	int e3=rectree->getNode(e1)->getFather()->getId();
	return rectree->getNode(e1)->getDist()-a1+dist(e3,0.0,e2,a2);
}

void ParamMR::makeDensity(ParamCons*p){
for (unsigned int k=0;k<edges.size();k++) {
while (p->getTree()->numRecEdge()>0) p->getTree()->remRecEdge(0);
for (unsigned int i=0;i<edges[k].size();i++) {
double afrom=edges[k][i].getTimeFrom(),ato=edges[k][i].getTimeTo();
int efrom=edges[k][i].getEdgeFrom(),eto=edges[k][i].getEdgeTo();
int start=edges[k][i].getStart(),end=edges[k][i].getEnd();
p->getTree()->addRecEdge(afrom,ato,start,end,efrom,eto);
}
p->account();
}
p->set();
}

RecEdge ParamMR::average(RecEdge*r1,RecEdge*r2,int coef)
{
    int    start=(coef*r1->getStart()+r2->getStart())/(coef+1);
    int    end  =(coef*r1->getEnd  ()+r2->getEnd  ())/(coef+1);
    int    efrom=0,eto=0;
    double afrom=0,ato=0;
    gotowards(r1->getEdgeFrom(),r1->getTimeFrom(),r2->getEdgeFrom(),r2->getTimeFrom(),&efrom,&afrom,dist(r1->getEdgeFrom(),r1->getTimeFrom(),r2->getEdgeFrom(),r2->getTimeFrom())/(coef+1.0));
    gotowards(r1->getEdgeTo  (),r1->getTimeTo  (),r2->getEdgeTo  (),r2->getTimeTo  (),&eto  ,&ato  ,dist(r1->getEdgeTo  (),r1->getTimeTo  (),r2->getEdgeTo  (),r2->getTimeTo  ())/(coef+1.0));
    RecEdge res(afrom,ato,start,end,efrom,eto);
    return res;
}

void ParamMR::gotowards(int e1,double a1,int e2,double a2,int*e3,double*a3,double d)
{
	if (e1==e2&&a1>=a2) {*e3=e1;*a3=a1-d;return;};
	if (e1==e2&&a1< a2) {*e3=e1;*a3=a1+d;return;};
	int cur=e2,cur2=0;
	while (cur!=-1&&cur!=e1) if (rectree->getNode(cur)->getFather()==NULL) cur=-1;else {cur2=cur;cur=rectree->getNode(cur)->getFather()->getId();}
	if (cur==-1) if (d<=rectree->getNode(e1)->getDist()-a1) {*e3=e1;*a3=a1+d;return;} else {gotowards(rectree->getNode(e1)->getFather()->getId(),0.0,e2,a2,e3,a3,d-rectree->getNode(e1)->getDist()+a1);return;};
	if (d<=a1) {*e3=e1;*a3=a1-d;return;} else {gotowards(cur2,rectree->getNode(cur2)->getDist(),e2,a2,e3,a3,d-a1);return;};
}

QString ParamMR::toString()
{
int n=data->getN()*2-1;
vector<vector<int> > s(n,vector<int>(n,0));
for (int i=0;i<rectree->numRecEdge();i++) {
  int efrom=rectree->getRecEdge(i)->getEdgeFrom();
  int eto  =rectree->getRecEdge(i)->getEdgeTo  ();
  if (rectree->getRecEdge(i)->getTimeFrom()==0) s[efrom][eto]=2; else s[efrom][eto]=1;
}
QString res;
for (int i=0;i<n;i++) for (int j=0;j<n;j++) res.append(QString("%1,").arg(s[i][j]));
res.chop(1);
return res;
}

void ParamMR::toCSV(QTextStream*out,int cutoff,int step)
{
  QString header;
  int n=getData()->getN()*2-1;
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) header.append(QString("\"from%1to%2\",").arg(i).arg(j));
  header.chop(1);
  *out<<header<<endl;
  int site=0;
  while (site<=l)
  {
  consensus(cutoff,site);
  *out<<toString()<<endl;
  site+=step;
  }
}

void ParamMR::toArtemis(QTextStream*out,int cutoff)
{
  consensusAllSites(cutoff);
  for (int i=0;i<rectree->numRecEdge();i++) {
    int    start=rectree->getRecEdge(i)->getStart()+1;
    int    end  =rectree->getRecEdge(i)->getEnd  ()+1;
    int    efrom=rectree->getRecEdge(i)->getEdgeFrom();
    int    eto  =rectree->getRecEdge(i)->getEdgeTo  ();
  *out<<"FT   misc_feature    "<<start*coef<<".."<<end*coef<<endl;
  if (rectree->getRecEdge(i)->getTimeFrom()==0) 
  *out<<"FT                   /note=\"recombination from "<<efrom<<" (or its children) "<<" to "<<eto<<"\""<<endl;
  else
  *out<<"FT                   /note=\"recombination from "<<efrom<<" to "<<eto<<"\""<<endl;
  }
  *out<<"ORIGIN"<<endl;
  QString str(l,QChar('N'));
  *out<<str<<endl;
}

void ParamMR::correctForPrior()
{
  double expected=getRho()*0.5*rectree->getTTotal()*getDelta()/(getDelta()*data->getB()+data->getL()-data->getB());
  int n=data->getN()*2-1;
  vector< vector<double> > corr(n,vector<double>(n,0.0));
  double s=0;
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
    double i0=rectree->getNode(i)->getAge();
    double j0=rectree->getNode(j)->getAge();
    double il=rectree->getNode(i)->getDist();
    double jl=rectree->getNode(j)->getDist();
    if (i==rectree->getN()*2-2) il=10.0;
    for (int a=0;a<100;a++) for (int b=0;b<100;b++)
      corr[i][j]+=rectree->priorEdge(i0+il*(a+1)/101.0,j0+jl*(b+1)/101.0);
    corr[i][j]*=jl*il/10000.0;
  s+=corr[i][j];
 }
 for (int i=0;i<n;i++) for (int j=0;j<n;j++) corr[i][j]*=expected/s;//Normalize
 for (int i=0;i<n;i++) for (int j=0;j<n;j++) for (int k=0;k<l;k++) states[i][j][k]-=corr[i][j]*its;
}
