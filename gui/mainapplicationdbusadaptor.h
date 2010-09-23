#include "mainwindowimpl.h"
#include <QDBusAbstractAdaptor>
 
 class MainApplicationAdaptor: public QDBusAbstractAdaptor
 {
     Q_OBJECT
     Q_CLASSINFO("weakarg", "org.gel.mauve.remote.WargInterface")

 private:
     MainWindowImpl *app;

 public:
     MainApplicationAdaptor(MainWindowImpl *window)
         : QDBusAbstractAdaptor(window), app(window)
     {
     }



 public slots:
     int getViewingSite()
     {
	return 0;
     }

     void setViewingSite(const int block, qlonglong site)
     {
	app->jumpToSite(site);
	return;
     }



 };
