#include "heatimpl.h"
//
HeatImpl::HeatImpl(int n, QWidget * parent, Qt::WFlags f)
    : QMainWindow(parent, f) {
  this->n=n;
  states=vector<vector<double> >(n,vector<double>(n,0.0));
  setupUi(this);
  its=0;
  expected=0.0;
}

void HeatImpl::account(ParamQt * p) 
{
  param=p;
  its++;
  RecTree * rectree=p->getTree();
  expected+=p->getRho()*0.5*rectree->getTTotal();
  if(!p->displaySet()) {
  for (int i=0;i<rectree->numRecEdge();i++)
    states[rectree->getRecEdge(i)->getEdgeFrom()][rectree->getRecEdge(i)->getEdgeTo()]++;
  }else{
	for (int i=0;i<rectree->numRecEdge();i++)
	{
	  vector<int> allto,allfrom;
	  allto =rectree->getAllSampledSeqs(rectree->getEdge(i)->getEdgeTo());
	  allfrom =rectree->getAllSampledSeqs(rectree->getEdge(i)->getEdgeFrom());
	  vector<int> allaffto=p->getDisplayTree()->getMinEdgeList(allto);
	  vector<int> allafffrom=p->getDisplayTree()->getMinEdgeList(allfrom);
	  for(unsigned int j=0;j<allaffto.size();j++) {
	    for(unsigned int k=0;k<allafffrom.size();k++) {
	      states[allafffrom[k]][allaffto[j]]++;
	    }
	  }
	}
  }
}

void HeatImpl::compute()
{
  table->setRowCount(n);
  table->setColumnCount(n);
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) 
    table->setItem(i,j,new QTableWidgetItem(QString::number(states[i][j]/its)));
}

void HeatImpl::compute_correct(int mode)
{
  RecTree * t=param->getTree();
  int num=0;for (int i=0;i<n;i++) for (int j=0;j<n;j++) num+=states[i][j];
  if (corr.size()==0) {
  corr=vector<vector<double> >(n,vector<double>(n,0.0));
  double s=0;
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
    double i0=t->getNode(i)->getAge();
    double j0=t->getNode(j)->getAge();
    double il=t->getNode(i)->getDist();
    double jl=t->getNode(j)->getDist();
    if (i==t->getN()*2-2) il=10.0;
    for (int a=0;a<100;a++) for (int b=0;b<100;b++)
      corr[i][j]+=t->priorEdge(i0+il*(a+1)/101.0,j0+jl*(b+1)/101.0);
    corr[i][j]*=jl*il/10000.0;
  s+=corr[i][j];
  }
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) corr[i][j]*=expected/s;
  }
  table->setRowCount(n);
  table->setColumnCount(n);
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) {
    QString val;
    if (mode==3) val=QString::number(corr[i][j]/its);else
    if (corr[i][j]==0.0||(corr[i][j]/its<3 && states[i][j]/its<3)) val="";else
    if (mode==2) val=QString::number(states[i][j]/corr[i][j]);
    else if (mode==1) val=QString::number((states[i][j]/its-corr[i][j]/its)/sqrt(corr[i][j]/its));
    table->setItem(i,j,new QTableWidgetItem(val));
  }
}

void HeatImpl::on_actionSave_as_activated() {
  if (actionGraph->isChecked()==false) {
      QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","CSV files (*.csv);;All files (*)");
      if (qstr==NULL) return;
      QFile file(qstr);
      if ( !file.open(QIODevice::WriteOnly)) return;
      QTextStream ts( &file );
      for (int i=0;i<table->rowCount();i++)
        for (int j=0;j<table->columnCount();j++) {
            if (table->item(i,j)!=NULL) ts<<table->item(i,j)->text();
            if (j<table->columnCount()-1) ts<<","; else ts<<endl;
          }
      file.close();
    } else {
      QString qstr = QFileDialog::getSaveFileName(this, tr("Save picture file"),".",tr("Joint Photographic Experts Group (*.jpg *.jpeg);;Windows Bitmap (*.bmp);;Portable Network Graphics (*.png);;Portable Pixmap (*.ppm);;X11 Bitmap (*.xbm *.xpm);;PostScript Format (*.ps);;Abode PDF Format (*.pdf)"));
      if (qstr==NULL) return;
      if (qstr.endsWith("ps") || qstr.endsWith("pdf")) {
          QPrinter qprint;
          qprint.setOutputFileName(qstr);
          qprint.setOrientation(QPrinter::Landscape);
          drawGraph(&qprint);
          return;
        }
      QImage image(width(),height(),QImage::Format_RGB32);
      image.invertPixels();//Fill image in white
      drawGraph(&image);
      image.save(qstr);
    }
    
}

void HeatImpl::on_actionQuit_activated() {
  close();
}

void HeatImpl::on_actionTable_activated() {
  actionTable->setChecked(true);
  actionGraph->setChecked(false);
  table->show();
  repaint();
}

void HeatImpl::on_actionGraph_activated() {
  actionGraph->setChecked(true);
  actionTable->setChecked(false);
  table->hide();
  repaint();
}

void HeatImpl::paintEvent(QPaintEvent*) {
  if (actionGraph->isChecked()) drawGraph(this);
}

void HeatImpl::drawGraph(QPaintDevice * qpd) {
  QPainter painter(qpd);
  double sup=0.0;
  for (int i=0;i<n;i++) for (int j=0;j<n;j++) sup=max(sup,abs(table->item(i,j)->text().toDouble()));
  if (!action0->isChecked()) sup=5.0;
  //Colorbar legend
  double leg=1.0;
  for (int i=0;i<11;i++) {
      painter.drawText(width()*0.92,height()*(0.11+0.8*i/10),QString::number(leg*sup,'f',3));
      leg-=0.2;
    }
  //YAxis labels
  for (int i=0;i<n;i+=5) {
      painter.drawText(0,height()*(0.1+0.8*i/n),width()*0.19,height()*0.8/n,Qt::AlignRight+Qt::AlignVCenter,QString::number(i));
      painter.drawText(width()*(0.2+0.6*i/n),0,width()*0.8/n,height()*0.17,Qt::AlignCenter+Qt::AlignTop,QString::number(i));
    }
  painter.translate(width()/5.0,height()/10.0);
  painter.scale(width()*0.6,height()*0.8);
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++) {
    	double l=table->item(i,j)->text().toDouble();
        int gain=255*l*(l>0)/sup;
        int loss=255*(-l)*(-l>0)/sup;
        if (gain>255) gain=255;
        if (loss>255) loss=255;
        painter.fillRect(QRectF(1.0*j/n,1.0*i/n,1.0/n,1.0/n),QColor(255-loss,255-gain-loss,255-gain));
      }
     }
  painter.drawRect(QRectF(0.0,0.0,1.0,1.0));
  for (int i=0;i<n;i+=5) {painter.drawLine(QLineF(1.0*i/n,0.0,1.0*i/n,1.0));
  painter.drawLine(QLineF(0.0,1.0*i/n,1.0,1.0*i/n));}
  
  //Color bar
  for (int i=-255;i<=255;i++) {
      int gain= i*(i>0);
      int loss=-i*(i<0);
      painter.fillRect(QRectF(1.1,1.0-1.0*(255+i)/511,0.05,1.0/511),QColor(255-loss,255-gain-loss,255-gain));
    }
  painter.drawRect(QRectF(1.1,0.0,0.05,1.0));
}

void HeatImpl::print(ostream* f_out)
{
  for (int i=0;i<n;i++) {
	for (int j=0;j<n-1;j++) *f_out<<table->item(i,j)->text().toStdString()<<",";
	 *f_out<<table->item(i,n-1)->text().toStdString()<<endl;
   }
}

void HeatImpl::on_menuPrior_correction_triggered(QAction * a)
{
  if (a==action0) {action0->setChecked(true );action1->setChecked(false);action2->setChecked(false);action3->setChecked(false);compute();}
  if (a==action1) {action0->setChecked(false);action1->setChecked(true );action2->setChecked(false);action3->setChecked(false);compute_correct(1);}
  if (a==action2) {action0->setChecked(false);action1->setChecked(false);action2->setChecked(true );action3->setChecked(false);compute_correct(2);}
  if (a==action3) {action0->setChecked(false);action1->setChecked(false);action2->setChecked(false);action3->setChecked(true );compute_correct(3);}
  repaint();
}
