#include "paramcons.h"
//
ParamCons::ParamCons(  ) 
{
	isCons=true;
	isSet=false;
	dp=NULL;
	dm=NULL;
	its=0;
}
//

ParamCons::~ParamCons()
{
}

void ParamCons::reset() 
{
	isSet=false;
	delete(dp);
	dp=NULL;
	delete(dm);
	dm=NULL;
	its=0;
}

void ParamCons::account(int id,bool getto)
{
	its++;
	if(!displayset|| displaytree==NULL) displaytree=rectree;
	if (dp==NULL) {dp=new DensityOnTree(displaytree,timeScale);}else dp->setTree(displaytree);
	if (dm==NULL) {dm=new DensityOnTree(displaytree,timeScale);}else dm->setTree(displaytree);
	for (int i=0;i<rectree->numRecEdge();i++)
	{
		if (gene>=0 && ((int)rectree->getRecEdge(i)->getStart()>=data->getBlocks()->at(gene) || (int)rectree->getRecEdge(i)->getEnd()<=data->getBlocks()->at(gene-1))) continue;
		if (id>=0 && getto && rectree->getRecEdge(i)->getEdgeTo()!=id) continue;
		else if(id>=0 && !getto && rectree->getRecEdge(i)->getEdgeFrom()!=id) continue;
		dp->add2(rectree,i,true);
		dm->add2(rectree,i,false);
	}
}

void ParamCons::accountrelative(int id,bool getto)
{
	its++;
	if(!displayset|| displaytree==NULL) displaytree=rectree;
	if (dp==NULL) dp=new DensityOnTree(displaytree,timeScale);else dp->setTree(displaytree);
	if (dm==NULL) dm=new DensityOnTree(displaytree,timeScale);else dm->setTree(displaytree);
	for (int i=0;i<rectree->numRecEdge();i++)
	{
		if (gene>=0 && ((int)rectree->getRecEdge(i)->getStart()>=data->getBlocks()->at(gene) || (int)rectree->getRecEdge(i)->getEnd()<=data->getBlocks()->at(gene-1))) continue;
		if (id>=0 && getto && rectree->getRecEdge(i)->getEdgeTo()!=id) continue;
		else if(id>=0 && !getto && rectree->getRecEdge(i)->getEdgeFrom()!=id) continue;

		dp->add2(rectree,i,true,rectree->getEdgeTreeTime(i));
		dm->add2(rectree,i,false,rectree->getEdgeTreeTime(i));
	}
}

void ParamCons::display(QPaintDevice * qpd) {
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
  
  if (dp!=NULL && dm!=NULL && isSet) {
  dp->display(&painter,x,y, rateScale/its,QBrush(Qt::red ));
  dm->display(&painter,x,y,-rateScale/its,QBrush(Qt::blue));}
  displayTree(&painter,x,y);
  delete(x);
  delete(y);
}
