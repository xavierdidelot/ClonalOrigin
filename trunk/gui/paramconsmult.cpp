#include "paramconsmult.h"
//
ParamConsMult::ParamConsMult(QStringList*nodes,QStringList*colors,bool denDep,bool colDep) 
{
	isCons=true;
	isSet=false;
	its=0;
	this->denDep=denDep;
	this->colDep=colDep;
	for (int i=0;i<nodes->size();i++) {this->nodes.push_back(nodes->at(i).toInt());
	brushes.push_back(QBrush(QColor(colors->at(i))));}
}
//

ParamConsMult::~ParamConsMult()
{
}

void ParamConsMult::account()
{
	its++;
	if (dps.size()==0) for (unsigned int i=0;i<nodes.size();i++) dps.push_back(new DensityOnTree(displaytree,timeScale));else for (unsigned int i=0;i<nodes.size();i++) dps[i]->setTree(displaytree);
	if (dms.size()==0) for (unsigned int i=0;i<nodes.size();i++) dms.push_back(new DensityOnTree(displaytree,timeScale));else for (unsigned int i=0;i<nodes.size();i++) dms[i]->setTree(displaytree);
	for (int i=0;i<rectree->numRecEdge();i++)
	{
		if (gene>=0 && ((int)rectree->getRecEdge(i)->getStart()>=data->getBlocks()->at(gene) || (int)rectree->getRecEdge(i)->getEnd()<=data->getBlocks()->at(gene-1))) continue;
		bool ok=false;
		for (unsigned int j=0;j<nodes.size();j++) {
		if (ok==false &&  colDep) ok=accept(rectree->getRecEdge(i)->getEdgeFrom(),nodes[j]);
		if (ok==false && !colDep) ok=accept(rectree->getRecEdge(i)->getEdgeTo  (),nodes[j]);
		if (ok==false) continue;
		dps[j]->add2(rectree,i,true);
		dms[j]->add2(rectree,i,false);
		}
	}
}

void ParamConsMult::display(QPaintDevice * qpd) {
  for (unsigned int i=6;i<dps.size();i++) brushes.push_back(brushes[i-6]);
  if(!displayset|| displaytree==NULL) displaytree=rectree;
  if (rectree==NULL) return;
  vector<double>*x,*y;
  getXY(&x,&y);
  QPainter painter(qpd);
  for (int i=0;i<displaytree->getN();i++)
    painter.drawText(QPointF(qpd->width()/10.0+x->at(i)*qpd->width()*0.7+3,qpd->height()/5.0+y->at(i)*qpd->height()*0.7+2),isolateName(i));
  if (timeScale<2.0) painter.drawText(qpd->width()*(0.45-0.05),qpd->height()*0.88,qpd->width()*0.1,qpd->height()*0.1,Qt::AlignCenter,QString::number(0.1));
  else painter.drawText(qpd->width()*(0.45-0.05),qpd->height()*0.875,qpd->width()*0.1,qpd->height()*0.1,Qt::AlignCenter,QString::number(1));
  double h=rateScale/displaytree->getN()/timeScale*0.005*displaytree->getNode(displaytree->getN()*2-2)->getAge();
  double mult=1.0;
  while (h*mult<0.0001) mult*=10.0;
  while (h*mult>0.001) mult/=10.0;
  painter.drawText(qpd->width()*0.11,qpd->height()*0.714,qpd->width()*0.1,qpd->height()*0.1,Qt::AlignCenter,QString::number(mult));
  painter.translate(qpd->width()/10.0,qpd->height()/5.0);
  painter.scale(qpd->width()*0.7,qpd->height()*0.7);
  double m=pow(h*mult,0.5);
  painter.fillRect(QRectF(0.05-m*0.5,0.8-m*0.5,m,m),QBrush(Qt::black));
  if (timeScale<2.0) {QLineF lin(0.5-0.05/timeScale,1.0,0.5+0.05/timeScale,1.0);painter.drawLine(lin);} else {QLineF lin(0.5-0.5/timeScale,1.0,0.5+0.5/timeScale,1.0);painter.drawLine(lin);}
  
  if (isSet) {
  for (int i=dps.size()-1;i>=0;i--) {
  if (!denDep) 
       dps[i]->display(&painter,x,y, -rateScale/its,brushes[i]);
  else dms[i]->display(&painter,x,y, -rateScale/its,brushes[i]);
  }}
  displayTree(&painter,x,y);
  delete(x);
  delete(y);
}
