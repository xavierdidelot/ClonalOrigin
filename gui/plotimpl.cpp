#include "plotimpl.h"
//
PlotImpl::PlotImpl( QWidget * parent, Qt::WFlags f)
    : QDialog(parent, f) {
  mode=0;
  setupUi(this);
  values=NULL;
  blocks=NULL;
}

PlotImpl::~PlotImpl() {
  delete(values);
}

void PlotImpl::paintEvent(QPaintEvent*) {
  if (mode==0) display_traj(this);
  if (mode==1) display_hist(this);
  if (mode==2) display_traj(this);
}

void PlotImpl::display_traj(QPaintDevice*qpd) {
  if (values==NULL) return;
  double minY=values->at(0);
  double maxY=values->at(0);
  if (blocks!=NULL) minY=0.0;
  for (unsigned int i=1;i<values->size();i++) {
      if (values->at(i)>maxY) maxY=values->at(i);
      if (values->at(i)<minY) minY=values->at(i);
    }
  double minX=0.0;
  double maxX=values->size()-1;
  QPainter painter(qpd);
  painter.drawText(QRectF(width()*0.0,height()*0.1,100,100),QString::number(maxY));
  painter.drawText(QRectF(width()*0.0,height()*0.9,100,100),QString::number(minY));
  painter.drawText(QRectF(width()*0.1,height()*0.95,100,100),QString::number(minX+1));
  painter.drawText(QRectF(width()*0.9,height()*0.95,100,100),QString::number(maxX+1));
  
  QPointF * points=(QPointF*)calloc(values->size(),sizeof(QPointF));
  for (unsigned int i=0;i<values->size();i++) {
      double x,y;
      x=i;
      y=values->at(i);
      x=(x-minX)/(maxX-minX);
      y=(y-minY)/(maxY-minY);
      points[i]=QPointF(x,1.0-y);
    }

  painter.translate(width()/10.0,height()/10.0);
  painter.scale(width()*0.8,height()*0.8);

  if (blocks!=NULL) {
    QPen pen(Qt::red);
    painter.setPen(pen);
    for (unsigned int i=1;i<blocks->size()-1;i++) {
      painter.fillRect(QRectF((blocks->at(i)-1-minX)/(maxX-minX),0.0,1.0/(maxX-minX),1.0),Qt::black);
      painter.drawLine(QLineF((blocks->at(i)-1-minX)/(maxX-minX),0.0,(blocks->at(i)-1-minX)/(maxX-minX),1.0));
     }
  }else meanAndCI();
    
  QPen pen(Qt::black);
  painter.setPen(pen);


  painter.drawPolyline(points,values->size());
  painter.drawRect(QRectF(-0.05,-0.05,1.1,1.1));
  free(points);
}

void PlotImpl::display_hist(QPaintDevice*qpd)
{
double minX=values->at(0);
double maxX=minX;
for (unsigned int i=1;i<values->size();i++) {
	if (values->at(i)<minX) minX=values->at(i);
	if (values->at(i)>maxX) maxX=values->at(i);
}
vector<int> hist=vector<int>(20,0);
for (unsigned int i=0;i<values->size();i++) hist[min(19.0,floor((values->at(i)-minX)/(maxX-minX)*20.0))]++;
int maxY=0;
for (unsigned int i=0;i<hist.size();i++) if (hist[i]>maxY) maxY=hist[i];
QPainter painter(qpd);
painter.drawText(QRectF(width()*0.0,height()*0.1,100,100),QString::number(maxY));
painter.drawText(QRectF(width()*0.0,height()*0.9,100,100),QString::number(0));
painter.drawText(QRectF(width()*0.1,height()*0.95,100,100),QString::number(minX));
painter.drawText(QRectF(width()*0.9,height()*0.95,100,100),QString::number(maxX));
painter.translate(width()/10.0,height()/10.0);
painter.scale(width()*0.8,height()*0.8);
for (unsigned int i=0;i<hist.size();i++)
  painter.drawRect(QRectF(1.0*i/20.0,1.0-1.0*hist[i]/maxY,1.0/20.0,1.0*hist[i]/maxY));
painter.drawRect(QRectF(-0.05,-0.05,1.1,1.1));
meanAndCI();
}

/*void PlotImpl::extractNumRecEdges(QDomDocument * domDoc) {
  QDomNodeList list=domDoc->elementsByTagName("Iteration");
  values=new vector<double>();
  int maxIteration=list.size();
  double mean=0.0;
  for (int i=0;i<maxIteration;i++) {      values->push_back(list.at(i).toElement().elementsByTagName("recedge").size());
      mean+=values->at(i);
    }
  mean/=values->size();
  vector<double> vals=vector<double>(*values);
  sort(vals.begin(),vals.end());
  double low= vals[(int)floor(0.025*vals.size())];
  double high=vals[(int)floor(0.975*vals.size())];
  setWindowTitle(QString::number(mean)+" ["+QString::number(low)+","+QString::number(high)+"]");
}

void PlotImpl::extractValues(QDomDocument * domDoc,string str) {
  QDomNodeList list=domDoc->elementsByTagName("Iteration");
  if (values!=NULL) delete(values);
  values=new vector<double>();
  int maxIteration=list.size();  
  for (int i=0;i<maxIteration;i++) 
    values->push_back(list.at(i).toElement().elementsByTagName(str.data()).at(0).firstChild().toText().data().toFloat());
  meanAndCI();
}*/

void PlotImpl::meanAndCI() {
double mean=0.0;
for (unsigned int i=0;i<values->size();i++) mean+=values->at(i);
mean/=values->size();
vector<double> vals=vector<double>(*values);
sort(vals.begin(),vals.end());
double low= vals[(int)floor(0.025*vals.size())];
double high=vals[(int)floor(0.975*vals.size())];
setWindowTitle(QString::number(mean)+" ["+QString::number(low)+","+QString::number(high)+"]");
}
/*
void PlotImpl::extractStat(QDomDocument * domDoc,string str,int b,int L) {
  if (QString::fromStdString(str).compare("rhotheta")==0) {
    extractValues(domDoc,"theta");
    vector<double> valTheta=vector<double>(*values);
    extractValues(domDoc,"rho");
    for (unsigned int i=0;i<values->size();i++) values->at(i)/=valTheta[i];
  } else {
    extractValues(domDoc,"delta");
    vector<double> valDelta=vector<double>(*values);
    extractValues(domDoc,"rho");
    for (unsigned int i=0;i<values->size();i++) values->at(i)*=valDelta[i]/(b*valDelta[i]+L-b);
    
  }
  meanAndCI();
}*/

void PlotImpl::on_exportButton_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save picture file"),".",tr("Joint Photographic Experts Group (*.jpg *.jpeg);;Windows Bitmap (*.bmp);;Portable Network Graphics (*.png);;Portable Pixmap (*.ppm);;X11 Bitmap (*.xbm *.xpm);;PostScript Format (*.ps);;Abode PDF Format (*.pdf)"));
  if (qstr==NULL) return;
    if (qstr.endsWith("ps") || qstr.endsWith("pdf")) {
  	QPrinter qprint;
	qprint.setOutputFileName(qstr);
	qprint.setOrientation( QPrinter::Landscape);
	  if (mode==0) display_traj(&qprint);
      if (mode==1) display_hist(&qprint);
      if (mode==2) display_traj(&qprint);
      return;
    }
  QImage image(width(),height(),QImage::Format_RGB32);
  image.invertPixels();//Fill image in white
  if (mode==0) display_traj(&image);
  if (mode==1) display_hist(&image);
  if (mode==2) display_traj(&image);
  image.save(qstr);
}
/*
void PlotImpl::extractGenome(QDomDocument * domDoc,vector<int> * bl,int which){
  QDomNodeList list=domDoc->elementsByTagName("Iteration");
  blocks=bl;
  values=new vector<double>(blocks->back(),0);
  for (int i=0;i<list.size();i++) {
    QDomNodeList listRec=list.at(i).toElement().elementsByTagName("recedge");
    for (int j=0;j<listRec.size();j++) {
    int start=listRec.at(j).toElement().elementsByTagName("start").at(0).firstChild().toText().data().toInt();
    int end  =listRec.at(j).toElement().elementsByTagName("end"  ).at(0).firstChild().toText().data().toInt();
    if (which==0) for (int k=start;k<end;k++) values->at(k)++;
    else if (which==1) values->at(start)++; 
    else values->at(end-1)++;
    }
   }
}*/

