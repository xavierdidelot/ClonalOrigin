#ifndef HEATIMPL_H
#define HEATIMPL_H
//
#include <QtGui>
#include "ui_heat.h"
#include <cmath>
#include "paramqt.h"
#include "../warg/src/param.h"

using namespace std;
//
class HeatImpl : public QMainWindow, public Ui::Heat
{
Q_OBJECT
public:
	HeatImpl(int n, QWidget * parent = 0, Qt::WFlags f = 0 );
	void account(ParamQt * p);
	void compute();
	void compute_correct(int mode);
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
        double expected;
        void paintEvent(QPaintEvent*);
        void drawGraph(QPaintDevice * qpd);
        vector<vector<double> > states;
        vector<vector<double> > corr;
};
#endif





