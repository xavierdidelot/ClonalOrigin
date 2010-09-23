#ifndef COLOUREDIMPL_H
#define COLOUREDIMPL_H
//
#include <QDialog>
#include "ui_coloured.h"
#include "mainwindowimpl.h"
//
class ColouredImpl : public QDialog, public Ui::Dialog
{
Q_OBJECT
public:
	ColouredImpl( QWidget * parent = 0, Qt::WFlags f = 0 );
private:
    void closeEvent(QCloseEvent * event);
private slots:
	void on_denArr_clicked();
	void on_colArr_clicked();
	void on_denDep_clicked();
	void on_colDep_clicked();
	void on_moveUp_clicked();
	void on_moveDown_clicked();
	void on_remove_clicked();
	void on_colour_clicked();
	void addGroup(int id);
	void on_buttonBox_accepted();
	void on_buttonBox_rejected();

};
#endif





