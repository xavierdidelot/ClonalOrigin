#ifndef PHEATIMPL_H
#define PHEATIMPL_H
//
#include <QtGui>
#include "ui_pheat.h"
#include <cmath>
#include "paramqt.h"
#include "../warg/src/param.h"

using namespace std;
//
class PHeatImpl : public QMainWindow, public Ui::PHeat
{
Q_OBJECT
public:
	PHeatImpl(int n, bool reldists=true, QWidget * parent = 0, Qt::WFlags f = 0 );
	void account(ParamQt * p);
	void compute(int mode);
	void print(ostream* f_out);
private slots:
	void on_menuPrior_correction_triggered(QAction * a);
	void on_actionSave_as_activated();
	void on_actionQuit_activated();
	void on_actionTable_activated();
	void on_actionGraph_activated();
protected:
        ParamQt * param;
	bool reldists;
        int n;
        int its;
        void paintEvent(QPaintEvent*);
        void drawGraph(QPaintDevice * qpd);
        vector<vector<int> > recmat;
        vector<vector<double> > priormat;
        vector<vector<double> > distmat;
};
#endif





