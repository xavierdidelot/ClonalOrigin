#ifndef PLOTIMPL_H
#define PLOTIMPL_H
//
#include "ui_plot.h"
#include <QtGui>
#include <QtXml>
#include <iostream>
#include <cmath>
#include <algorithm>
//

using namespace std;

class PlotImpl : public QDialog, public Ui::Plot {
      Q_OBJECT
    public:
      PlotImpl( QWidget * parent = 0, Qt::WFlags f = 0 );
      ~PlotImpl();
      //void extractValues(QDomDocument * domDoc,string str);
      //void extractStat(QDomDocument * domDoc,string str,int b=1,int L=1);
      void meanAndCI();
      //void extractNumRecEdges(QDomDocument * domDoc);
      //void extractGenome(QDomDocument * domDoc,vector<int> * bl,int which);
      inline void setValues(vector<double>*v) {values=v;}
      void trajectory();
      void hist();
      inline void setMode(int m) {mode=m;}
      inline void setBlocks(vector<int>*bl) {blocks=bl;}
    private slots:
      void on_exportButton_clicked();
    protected:
      void paintEvent(QPaintEvent*);
      void display_traj(QPaintDevice*qpd);
      void display_hist(QPaintDevice*qpd);
      vector<double> * values;
      vector<int> * blocks;
      int mode;
  };
#endif
