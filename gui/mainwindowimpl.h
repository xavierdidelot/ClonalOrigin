#ifndef MAINWINDOWIMPL_H
#define MAINWINDOWIMPL_H
//
#include "ui_mainwindow.h"
#include "plotimpl.h"
#include "gelmanrubinimpl.h"
#include "paramcons.h"
#include "paramconsmult.h"
#include <QtSvg>
#include "outputfile.h"
#include "parammr.h"
#include "paramtreecons.h"
#include "colouredimpl.h"
#include "heatimpl.h"
#include "pdimpl.h"
#include "pheatimpl.h"
//
class MainWindowImpl : public QMainWindow, public Ui::MainWindow {
      Q_OBJECT
    public:
      MainWindowImpl( QWidget * parent = 0, Qt::WFlags f = 0 );
      virtual ~MainWindowImpl();
      void openXMLFile(QStringList qstrs);
      void doColourPlot(QStringList*nodes,QStringList*colors,bool denDep,bool colDep);
	void on_actionExport_in_CSV_format_activated(QString qstr,int grmode=2);
	void on_actionScore_against_true_tree_activated(QString qstr,QString dest);
	void on_actionTest_convergence_of_trees_activated(QStringList files,QString dest);
	void on_actionHeat_map_activated(int correctforprior,QString dest);
	void on_actionPd_map_activated(int correctforprior,QString dest);
	void on_actionPheat_map_activated(int correctforprior,QString dest,bool reldists=true);
	void on_actionGelman_Rubin_test_activated(QStringList files,QString dest,int grmode=2);
	void on_actionMajority_rule_consensus_of_trees_activated(int p,QString dest);
	void on_actionExtended_consensus_of_trees_activated(QString dest);
        void on_actionCombine_output_files_activated(QStringList files,QString dest);
	void on_actionThinIterations_activated(QString outputfile,int thin);
	void on_actionExtractIterations_activated(QString outputfile,int istart,int iend);
	void jumpToSite(int site);
    inline void setRecView(int i){param->setRecView(i);}	

    signals:
      void wasDbClicked(int node,QMouseEvent * event);
    private slots:
      void computeAfterDbClick(int id,QMouseEvent * event);
      void on_actionGelman_Rubin_test_activated();
      void on_actionShowComment_activated();
      void on_actionExit_activated();
      void on_actionSave_picture_activated();
      void on_actionReOpen_output_file_activated();
      void on_actionOpen_output_file_activated();
      void on_actionAbout_activated();
      void on_menuPlot_triggered(QAction* action);
      void on_menuVisualisation_triggered(QAction* action);
      void on_actionCombine_output_files_activated();
      void on_actionExtractIterations_activated();
      void on_actionThinIterations_activated();
	void on_actionCompute_prior_corrected_activated();
	void on_actionExport_to_Artemis_activated();
	void on_actionExport_CSV_activated();
	void on_actionExport_movie_activated();
	void on_actionJump_to_site_activated();
	void on_actionSet_cutoff_activated();
	void on_actionPrev_site_activated();
	void on_actionNext_site_activated();
	void on_actionAll_genes_activated();
	void on_actionOnly_one_gene_activated();
	void on_actionSave_movie_activated();
	void on_actionToggleView_activated();
	void on_actionShowCF_activated();
	void on_actionDeparture_arrival_density_activated();
	void on_actionMajority_rule_Consensus_activated();
	void on_actionExport_as_Tracer_file_activated();
	void on_actionSave_tree_to_Newick_file_activated();
	void on_actionSave_trees_to_Nexus_file_activated();
	void on_actionMajority_rule_Density_activated();
	void on_actionRelative_Departure_Arrival_density_activated();
	void on_actionColour_density_plot_activated();
	void on_actionHeat_map_activated();
	void on_actionPd_map_activated();
	void on_actionPheat_map_activated();
	void on_actionPheat_map_rel_activated();
	void on_actionNameFormatm1_activated();
	void on_actionNameFormat0_activated();
	void on_actionNameFormat1_activated();
	void on_actionNameFormat2_activated();
	inline void on_actionMajority_rule_consensus_of_trees_activated(){
		on_actionMajority_rule_consensus_of_trees_activated(-1,QString());
	}
	inline void on_actionExtended_consensus_of_trees_activated(){
		on_actionExtended_consensus_of_trees_activated(QString());
	}
	inline void on_actionExport_in_CSV_format_activated(){
		on_actionExport_in_CSV_format_activated(QString());
	}
	inline void on_actionScore_against_true_tree_activated(){
		on_actionScore_against_true_tree_activated(QString(),QString());
	}
	inline void on_actionTest_convergence_of_trees_activated(){
		on_actionTest_convergence_of_trees_activated(QStringList(),QString());
	};
    int doDeparture_arrival_density(ParamCons * paramcons,bool relative=false);
    ParamTreeCons * newExtConsTree();///<returns a pointer to a new extended consensus tree
    void computeExplorer(bool correct);
	inline void checkDataLoaded(){
	  	if(data==NULL) {QStringList qstrs = QFileDialog::getOpenFileNames(this, tr("Open data file"),".","XMFA files (*.xmfa);;All files (*)");
		data=new Data(qstrs[0].toStdString()); // we have to keep reloading the data
		WargXml infile(outputfilenames[0].toStdString());
		data->subset(outputFile->getRegions(),-1);
		param->setData(data);
		}
	}
    protected:
      QStringList outputfilenames;
      ParamQt * param;
      ParamQt * ioparam;
      int explorerSite;
      int explorerCutoff;
      //QDomDocument * domDoc;
      //QDomNodeList its;
      OutputFile * outputFile;
      Data * data;
      //int currentIteration;
      //int maxIteration;
      void paintEvent(QPaintEvent*);
      void loadIteration(bool startOver=false);
      void mouseDoubleClickEvent ( QMouseEvent * event );
      void displayStatus(QString str=tr(""));
  };
#endif
