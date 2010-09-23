#include "colouredimpl.h"
//
ColouredImpl::ColouredImpl( QWidget * parent, Qt::WFlags f) 
	: QDialog(parent, f)
{
	setAttribute(Qt::WA_DeleteOnClose);
	setupUi(this);
}
//

void ColouredImpl::on_buttonBox_accepted()
{
	QStringList * nodes=new QStringList();
	QStringList * colors=new QStringList();
	for (int i=0;i<listWidget->count();i++) {
		QListWidgetItem * item=listWidget->item(i);
		nodes->append(item->text());
		colors->append(item->textColor().name());
	}
	if (nodes->size()>0) ((MainWindowImpl*)parent())->doColourPlot(nodes,colors,denDep->isChecked(),colDep->isChecked());
    close();
}

void ColouredImpl::on_buttonBox_rejected()
{
    close();
}

void ColouredImpl::closeEvent(QCloseEvent * event)
{
	event=0;
	connect(((MainWindowImpl*)parent()),SIGNAL(wasDbClicked(int,QMouseEvent*)),((MainWindowImpl*)parent()),SLOT(computeAfterDbClick(int,QMouseEvent*)));
    destroy();
}

void ColouredImpl::addGroup(int id)
{
	QListWidgetItem * item=new QListWidgetItem(QString::number(id));
	item->setTextColor(Qt::red);
	listWidget->addItem(item);
}

void ColouredImpl::on_colour_clicked()
{
	QListWidgetItem * item=listWidget->currentItem();
	if (item==NULL) return;
	QColor c0=item->textColor();
	QColor c=QColorDialog::getColor(c0);
	if (!c.isValid()) return;
	item->setTextColor(c);
}

void ColouredImpl::on_moveUp_clicked()
{
	QListWidgetItem * item=listWidget->currentItem();
	if (item==NULL) return;
	int c=listWidget->currentRow();
	if (c==0) return;
	item=listWidget->takeItem(listWidget->currentRow());
	listWidget->insertItem(c-1,item);
	listWidget->setCurrentRow(c-1);
}

void ColouredImpl::on_moveDown_clicked()
{
	QListWidgetItem * item=listWidget->currentItem();
	if (item==NULL) return;
	int c=listWidget->currentRow();
	if (c==listWidget->count()-1) return;
	item=listWidget->takeItem(listWidget->currentRow());
	listWidget->insertItem(c+1,item);
	listWidget->setCurrentRow(c+1);
}

void ColouredImpl::on_remove_clicked()
{
	QListWidgetItem * item=listWidget->currentItem();
	if (item==NULL) return;
	listWidget->takeItem(listWidget->currentRow());
}

void ColouredImpl::on_denArr_clicked()
{
	denDep->setChecked(false);
}

void ColouredImpl::on_colArr_clicked()
{
	colDep->setChecked(false);
}

void ColouredImpl::on_denDep_clicked()
{
	denArr->setChecked(false);
}

void ColouredImpl::on_colDep_clicked()
{
	colArr->setChecked(false);
}
