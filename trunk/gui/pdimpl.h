#ifndef PDIMPL_H
#define PDIMPL_H
//
#include <QtGui>
#include "ui_pd.h"
#include <cmath>
#include "paramqt.h"
#include "../warg/src/param.h"

using namespace std;
//
class PdImpl : public QMainWindow, public Ui::Pd
{
Q_OBJECT
public:
	PdImpl(int n, int dist=1, QWidget * parent = 0, Qt::WFlags f = 0 );
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
        int n;
        int its;
	int sitedist;
        void paintEvent(QPaintEvent*);
        void drawGraph(QPaintDevice * qpd);
        vector<vector<double> > recmat;
        vector<vector<double> > treemat;
};
#endif





