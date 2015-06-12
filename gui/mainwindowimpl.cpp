#include "mainwindowimpl.h"
#include <QDBusConnection>
#include <QDBusMessage>

//
MainWindowImpl::MainWindowImpl( QWidget * parent, Qt::WFlags f):QMainWindow(parent, f) {
  param=new ParamQt();
  outputFile=NULL;
  data=NULL;
  explorerCutoff=70;
  connect(this,SIGNAL(wasDbClicked(int,QMouseEvent*)),this,SLOT(computeAfterDbClick(int,QMouseEvent*)));
  setupUi(this);
}

MainWindowImpl::~MainWindowImpl() {
  if(param) delete(param);
  if(outputFile) delete(outputFile);
}

void MainWindowImpl::on_actionAbout_activated() {
  QMessageBox::about(0, "About...","This is weakarg-gui version 0.1.");
}

void MainWindowImpl::paintEvent(QPaintEvent*) {
  param->display(this);
}

void MainWindowImpl::on_actionReOpen_output_file_activated() {
  outputFile->startOver();
  outputFile->reset();
  param->setNames(outputFile->getNames());
  param->setBlocks(outputFile->getBlocks());
  loadIteration();
  repaint();
}

void MainWindowImpl::on_actionOpen_output_file_activated() {
  QStringList qstrs = QFileDialog::getOpenFileNames(this, tr("Open output file(s)"),".","XML files (*.xml);;All files (*)");
  if (qstrs.isEmpty()) return;
  openXMLFile(qstrs);
}

void MainWindowImpl::openXMLFile(QStringList qstrs) {
  outputfilenames=qstrs;
  if(outputFile) delete(outputFile);
  outputFile=new OutputFile(qstrs,false);//set this to true to recover the previous behaviour of initialising on loading.
  param->setBlocks(outputFile->getBlocks());
  param->setNames(outputFile->getNames());
  param->setLabels(QStringList());
  param->clearTreeData();
  data=NULL;
  //if(data!=NULL) delete(data);
  loadIteration();
  repaint();
}

void MainWindowImpl::loadIteration(bool startOver) {
  if (outputFile==NULL) return;
  if (param->isCons) {ParamQt * param2=new ParamQt();
  param2->setTree(param->getTree());
  param2->setBlocks(outputFile->getBlocks());
  param2->setNames(outputFile->getNames());
  if(param) delete(param);
  param=param2;}
  if (startOver) {outputFile->startOver();  param->setNames(outputFile->getNames());}
  outputFile->getIt(param);
  displayStatus();
}

void MainWindowImpl::on_actionExit_activated() {
  close();
}

void MainWindowImpl::on_actionShowComment_activated() {
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  QMessageBox::about(0,"Outputfile comment",outputFile->getComment());
}

void MainWindowImpl::on_actionSave_picture_activated() {
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save picture file"),".",tr("Joint Photographic Experts Group (*.jpg *.jpeg);;Windows Bitmap (*.bmp);;Portable Network Graphics (*.png);;Portable Pixmap (*.ppm);;X11 Bitmap (*.xbm *.xpm);;SVG Format (*.svg);;PostScript Format (*.ps);;Abode PDF Format (*.pdf)"));
  if (qstr==NULL) return;
  if (qstr.endsWith("svg")) {
  	QSvgGenerator qsvg;
	qsvg.setFileName(qstr);
	qsvg.setSize(QSize(width(),height()));
      param->display(&qsvg);
      return;
    }
    if (qstr.endsWith("ps") || qstr.endsWith("pdf")) {
  	QPrinter qprint;
	qprint.setOutputFileName(qstr);
//	qprint.setOrientation( QPrinter::Landscape);
      param->display(&qprint);
      return;
    }
  QImage image(width(),height(),QImage::Format_ARGB32);
  image.invertPixels();//Fill image in white
  param->display(&image);
  image.save(qstr);
}

void MainWindowImpl::on_actionSave_trees_to_Nexus_file_activated() {
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save trees file"),".","Nexus file (*.nex);;All files (*)");
  if (qstr==NULL) return;
  QFile file(qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream ts( &file );
  ts<<"#nexus"<<endl<<"begin trees;"<<endl<<"  translate"<<endl;
  for(int i=0;i<param->getTree()->getN()-1;i++) ts<<"    "<<i<<" "<<i<<","<<endl;
  ts<<"    "<<param->getTree()->getN()-1<<" "<<param->getTree()->getN()-1<<endl<<"  ;"<<endl;
  outputFile->startOver();
  int treeon=0;
  while (outputFile->getIt(param))
  {
    displayStatus();
    ts << "  tree "<<treeon++<<" = "<<"  "<<QString::fromStdString(param->getTree()->newickNoInternalLabels())<<endl;
  }
  ts<<"end;"<<endl;
  file.close();



}


void MainWindowImpl::on_actionSave_tree_to_Newick_file_activated() {
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save tree file"),".","Newick file (*.nwk *.newick *.tree);;All files (*)");
  if (qstr==NULL) return;
  QFile file(qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream ts( &file );
  ts << QString::fromStdString(param->getTree()->newick());
  file.close();
}

void MainWindowImpl::on_actionGelman_Rubin_test_activated(QStringList files,QString dest,int grmode) {
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  if(files.size()==0) files  = QFileDialog::getOpenFileNames(this, tr("Select File(s)"),".","XML files (*.xml);;All files (*)");
  if (files.size()==0) return;
  GelmanRubinImpl * gr=new GelmanRubinImpl(this);
  if(dest.length()==0) {
    gr->showOptions();
    gr->setFiles(param,outputFile,files);
//    gr->compute(param,outputFile,&files,NULL,true,true,true);
    gr->show();
  }else {//gr->compute(param,outputFile,&files,&cout,true,true,true);
  //}else{
	ostream* f_out;
	if(dest.compare("COUT")==0) f_out=&cout;
	else f_out = new ofstream(qPrintable(dest));
	if(grmode==1) gr->compute(param,outputFile,&files,f_out,true,false,false);
	else if(grmode==2) gr->compute(param,outputFile,&files,f_out,true,true,false);
	else if(grmode==3) gr->compute(param,outputFile,&files,f_out,true,false,true);
	else gr->compute(param,outputFile,&files,f_out,true,true,true);
	if(dest.compare("COUT")!=0) delete f_out;
  }
  loadIteration();
}

void MainWindowImpl::on_actionGelman_Rubin_test_activated(){ 
on_actionGelman_Rubin_test_activated(QStringList(),QString(""));
}


void MainWindowImpl::on_actionExport_as_Tracer_file_activated()
{
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","Tracer files (*.log);;All files (*)");
  if (qstr==NULL) return;
  GelmanRubinImpl * gr=new GelmanRubinImpl(this);
  QStringList strlist=QStringList(qstr);
    gr->setFiles(param,outputFile,strlist);
    gr->setExport();
    gr->setSep(string("\t"));
    gr->show();
 // gr->outputTracer(param,outputFile,&qstr);
  outputFile->startOver();
  loadIteration();
  //delete(gr);
}

void MainWindowImpl::on_actionExport_in_CSV_format_activated(QString qstr,int grmode)
{
  bool vis=false;
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  if(qstr.length()==0) {vis=true;qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","CSV index files (*.csv);;All files (*)");}
  if (qstr==NULL) return;
  GelmanRubinImpl * gr=new GelmanRubinImpl(this);
  if(vis){
  QStringList strlist=QStringList(qstr);
    gr->setFiles(param,outputFile,strlist);
    gr->setExport();
    gr->show();
  }else {
    if(grmode==1) gr->outputTracer(param,outputFile,&qstr,true,true,false,false);
    else if(grmode==2) gr->outputTracer(param,outputFile,&qstr,true,true,true,false);
    else if(grmode==3) gr->outputTracer(param,outputFile,&qstr,true,true,false,true);
    else gr->outputTracer(param,outputFile,&qstr,true,true,true,true);
    delete(gr);
  }
  outputFile->startOver();
  loadIteration();
}

void MainWindowImpl::on_menuPlot_triggered(QAction* a) {
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
// Should give the user a warning this is slow!
  if(!outputFile->isinitialised()) outputFile->countVectors();
  PlotImpl * pi=new PlotImpl(this);
  if (a==action11 || a==action12) pi->setValues(outputFile->getThetas());
  if (a==action21 || a==action22) pi->setValues(outputFile->getRhos());
  if (a==action31 || a==action32) pi->setValues(outputFile->getDeltas());
  if (a==action41 || a==action42) pi->setValues(outputFile->getNumRecEdges());
  if (a==action51 || a==action52) pi->setValues(outputFile->getLikelihoods());
  if (a==action61 || a==action62) pi->setValues(outputFile->getRhoOverTheta());
  if (a==action71 || a==action72) {
	checkDataLoaded();
	pi->setValues(outputFile->getRoverM(param));
	repaint();
  }
  if (a==action81 || a==action82) pi->setValues(outputFile->getTMRCA());
  if (a==action91 || a==action92) pi->setValues(outputFile->getTTotal());
  if (a==actionA1 || a==actionA2) pi->setValues(outputFile->getPriors());
  if (a==actionB1 || a==actionB2) pi->setValues(outputFile->getPosteriors());
  if (a==actionC1 || a==actionC2) pi->setValues(outputFile->getThetaPerSite());
  if (a==actionD1 || a==actionD2) pi->setValues(outputFile->getRhoPerSite());
  if (a==actionRecombination) pi->setValues(outputFile->getGenoRec());
  if (a==actionRelRecombination) pi->setValues(outputFile->getRelGenoRec(param));
  if (a==actionStarting_points) pi->setValues(outputFile->getGenoBeg());
  if (a==actionEnding_points) pi->setValues(outputFile->getGenoEnd());
  if (a==action11 || a==action21 || a==action31 || a==action41 || a==action51 || a==action61 || a==action71 || a==action81 || a==action91 || a==actionA1 || a==actionB1 || a==actionC1 || a==actionD1) pi->setMode(0);else 
  if (a==action12 || a==action22 || a==action32 || a==action42 || a==action52 || a==action62 || a==action72 || a==action82 || a==action92 || a==actionA2 || a==actionB2 || a==actionC2 || a==actionD2) pi->setMode(1);else
  {pi->setBlocks(param->getData()->getBlocks());pi->setMode(2);}
  pi->show();
}

void MainWindowImpl::on_menuVisualisation_triggered(QAction* a) {
  if (a==actionNext_iteration) loadIteration();
  if (a==actionFirst_iteration) loadIteration(true);
  if (a==actionChange_scale) {
  bool ok;
  double r=QInputDialog::getDouble(this,"Enter rate scale","Enter value for the scale of the rates:",param->getRateScale(),0,2147483647,10,&ok);
  if (!ok) return;
    param->setRateScale(r);
    displayStatus(tr("Scale=")+QString::number(param->getRateScale()));
  };
  if (a==actionTime_scale_up) param->incrTimeScale();
  if (a==actionTime_scale_down) param->decrTimeScale();
  //if (a==actionTime_scale_up || a==actionTime_scale_down)
  displayStatus();
  repaint();
}

void MainWindowImpl::on_actionCombine_output_files_activated(QStringList files,QString dest)
{
  QFile file(dest);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream out(&file);
  out << "<?xml version = '1.0' encoding = 'UTF-8'?>" << endl;
  out << "<outputFile>" << endl;
  for (int i=0;i<files.size();i++) {
  QFile file2(files[i]);
  if ( !file2.open(QIODevice::ReadOnly)) return;
  while (!file2.atEnd()) {
      QByteArray line = file2.readLine();
      if (!line.contains("xml version") && !line.contains("outputFile")) out<<line;
  }
  file2.close();
  }
  out<<"</outputFile>"<<endl;
  file.close();
  //Open the result
  openXMLFile(QStringList(dest));
}

void MainWindowImpl::on_actionCombine_output_files_activated()
{
  //Select files to combine
  QStringList files  = QFileDialog::getOpenFileNames(this, tr("Select File(s)"),".","XML files (*.xml);;All files (*)");
  if (files.size()==0) return;
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","XML files (*.xml);;All files (*)");
  if (qstr==NULL) return;
  //Combine them
  on_actionCombine_output_files_activated(files,qstr);
 
}

void MainWindowImpl::on_actionExtractIterations_activated(QString outputfile,int istart,int iend)
{
  if (outputfile==NULL) return;
  //Extract them
  QFile file(outputfile);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream out(&file);
  QFile filefrom(outputFile->getFileName());
  if ( !filefrom.open(QIODevice::ReadOnly)) return;
  int iton=-1;
  while (!filefrom.atEnd()) {
      QByteArray line = filefrom.readLine();
      if(line.contains("<Iteration>")) iton++;
      if (iton<0 || line.contains("outputFile") || (iton>=istart && iton<=iend)) out<<line;
  }
  filefrom.close();
  file.close();
  //Open the result
  openXMLFile(QStringList(outputfile));
}

void MainWindowImpl::on_actionExtractIterations_activated()
{
  //Select iterations to combine
	 bool ok;
  int istart = QInputDialog::getInteger(this, tr("Start Sample Number"),tr("Enter the mininum sample number:"), 0, 0, 2147483647, 1, &ok);
  int iend = QInputDialog::getInteger(this, tr("End Sample Number"),tr("Enter the maximum sample number:"), 0, 0, 2147483647, 1, &ok);
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","XML files (*.xml);;All files (*)");
  on_actionExtractIterations_activated(qstr,istart,iend);
}

void MainWindowImpl::on_actionThinIterations_activated(QString outputfile,int thin)
{
  if (outputfile==NULL) return;
  //Extract them
  QFile file(outputfile);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream out(&file);
  QFile filefrom(outputFile->getFileName());
  if ( !filefrom.open(QIODevice::ReadOnly)) return;
  int iton=-1;
  while (!filefrom.atEnd()) {
      QByteArray line = filefrom.readLine();
      if(line.contains("<Iteration>")) iton++;
      if (iton<0 || line.contains("outputFile") || (iton%thin==0)) out<<line;
  }
  filefrom.close();
  file.close();
  //Open the result
  openXMLFile(QStringList(outputfile));
}

void MainWindowImpl::on_actionThinIterations_activated()
{
  //Select iterations to combine
	 bool ok;
  int thin = QInputDialog::getInteger(this, tr("Thin Amount"),tr("Enter the thinning step:"), 0, 0, 2147483647, 1, &ok);
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","XML files (*.xml);;All files (*)");
  on_actionThinIterations_activated(qstr,thin);
}

void MainWindowImpl::on_actionAll_genes_activated()
{
	actionAll_genes->setChecked(true);
	actionOnly_one_gene->setChecked(false);
	param->setGene(-1);
	repaint();
}

void MainWindowImpl::on_actionOnly_one_gene_activated()
{
	actionAll_genes->setChecked(false);
	actionOnly_one_gene->setChecked(true);
	 bool ok;
     int i = QInputDialog::getInteger(this, tr("Pick a gene"),tr("Show events for gene number:"), 1, 1, param->getData()->getB(), 1, &ok);
     if (ok) param->setGene(i);
     repaint();
}

void MainWindowImpl::on_actionToggleView_activated(){
	param->toggleRecView();
	if(param->getRecView()==0)statusBar()->showMessage("Showing all recombination equally.");
	else if(param->getRecView()==1)statusBar()->showMessage("Showing recombination weighted by tract length.");
	else if(param->getRecView()==1)statusBar()->showMessage("Not showing recombination.");
	repaint();
}

void MainWindowImpl::on_actionShowCF_activated(){
	checkDataLoaded();
	outputFile->makeCF(param);
	param->setRecView(2);
	repaint();	
}

void MainWindowImpl::on_actionSave_movie_activated()
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save movie file"),".",tr("WMV files (*.wmv);;AVI files (*.avi)"));
  if (qstr==NULL) return;
  QImage image(width(),height(),QImage::Format_ARGB32);
  outputFile->startOver();
  QString time=QTime::currentTime().toString();
  while (outputFile->getIt(param))
  {
  image.fill(0);
  image.invertPixels();//Fill image in white
  param->display(&image);
  image.save(QString("/tmp/%1it%2.jpg").arg(time).arg(QString::number(outputFile->getCurIt()),7,'0'));
  }
  QProcess ffmpeg;
  if (qstr.endsWith("wmv")) 
  ffmpeg.start("ffmpeg", QString("-y -r 10 -i /tmp/%1it%07d.jpg -vcodec wmv1 %2").arg(time).arg(qstr).split(" "));
  else 
  ffmpeg.start("ffmpeg", QString("-y -r 10 -i /tmp/%1it%07d.jpg %2").arg(time).arg(qstr).split(" "));
  ffmpeg.waitForFinished();
}

void MainWindowImpl::on_actionDeparture_arrival_density_activated()
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  ParamCons * paramcons=new ParamCons();
  paramcons->setBlocks(outputFile->getBlocks());
  paramcons->setNames(outputFile->getNames());
  paramcons->setTimeScale(param->getTimeScale());
  paramcons->setGene(param->getGene());
  if(param->displaySet()) paramcons->newDisplayTree(param->getDisplayTree(),false);
  if(param) delete(param);
  param=paramcons;
  if(doDeparture_arrival_density(paramcons)==-1) doDeparture_arrival_density(paramcons);
  paramcons->set();
  repaint();
}

int MainWindowImpl::doDeparture_arrival_density(ParamCons * paramcons,bool relative)
{
  outputFile->startOver();
  bool firstit=true;
  while (outputFile->getIt(paramcons))
    {
     if(firstit && !paramcons->displaySet()) {
	  paramcons->newDisplayTree(paramcons->getTree());
	  firstit=false;
      }
      if(!paramcons->displaySet()) {
	string s1(paramcons->getDisplayTree()->newick());
	string s2(paramcons->getTree()->newick());
	if(s1.compare(s2)!=0) {
		ParamTreeCons * paramtreecons=newExtConsTree();
		paramcons->unsetDisplayTree();
		paramcons->newDisplayTree(paramtreecons->getTree());
		delete(paramtreecons);
		repaint();
		return -1;
	}
      }
      displayStatus();
      if(relative)paramcons->accountrelative();
      else paramcons->account();
    }
    return(0);
}


void MainWindowImpl::on_actionRelative_Departure_Arrival_density_activated()
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  ParamCons * paramcons=new ParamCons();
  paramcons->setBlocks(outputFile->getBlocks());
  paramcons->setNames(outputFile->getNames());
  paramcons->setTimeScale(param->getTimeScale());
  paramcons->setGene(param->getGene());
  if(param) delete(param);
  param=paramcons;

  if(doDeparture_arrival_density(paramcons)==-1) doDeparture_arrival_density(paramcons,true);

  paramcons->set();
  repaint();
}

void MainWindowImpl::on_actionMajority_rule_Consensus_activated()
{
  computeExplorer(false);
}
 
void MainWindowImpl::computeExplorer(bool correct)
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  // tell Mauve to order its sequences and focus on site 0
  // FIXME: need to figure out how to get the genome indices and send them as an array.
  QDBusMessage m = QDBusMessage::createMethodCall("org.gel.mauve.remote.MauveInterface","/MauveInterface","","hackOrder");
  QDBusMessage response = QDBusConnection::sessionBus().call(m);
  m = QDBusMessage::createMethodCall("org.gel.mauve.remote.MauveInterface","/MauveInterface","","setDisplayBlockAndColumn");
  m << (int)0;
  m << (qlonglong)0;
  m << (qlonglong)500;
  response = QDBusConnection::sessionBus().call(m);

  explorerSite=0;
  ParamMR * paramMR=new ParamMR();
  paramMR->setBlocks(outputFile->getBlocks());
  paramMR->setNames(outputFile->getNames());
  paramMR->setTimeScale(param->getTimeScale());
  if(param) delete(param);
  param=paramMR;
  outputFile->startOver();
  bool firstit=true;
  while (outputFile->getIt(paramMR))
  {
      if(firstit) {
	  paramMR->newDisplayTree(paramMR->getTree());
	  firstit=false;
      }
    //string s1(paramMR->getDisplayTree()->newick());
	//string s2(paramMR->getTree()->newick());
	//if(s1.compare(s2)!=0) {QMessageBox::about(0,"Information","Cannot do majority rule for non-constant trees.");return;}
    displayStatus();
    paramMR->account();
  }
  if (correct) paramMR->correctForPrior();
  paramMR->consensus(explorerCutoff);
  repaint();
}

void MainWindowImpl::on_actionMajority_rule_Density_activated()
{
  on_actionMajority_rule_Consensus_activated();
  ParamCons * paramcons=new ParamCons();
  paramcons->setBlocks(outputFile->getBlocks());
  paramcons->setNames(outputFile->getNames());
  paramcons->setTimeScale(param->getTimeScale());
  outputFile->startOver();
  outputFile->getIt(paramcons);
  ((ParamMR*)param)->makeDensity(paramcons);
  if(param) delete(param);
  param=paramcons;
  repaint();
}


void MainWindowImpl::on_actionMajority_rule_consensus_of_trees_activated(int cutoff, QString dest)
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  bool ok=true;
  if(cutoff<0) cutoff=QInputDialog::getInteger(this,"Enter cutoff","Enter value of cutoff for majority-rule consensus:",95,0,100,1,&ok);
  if (!ok || cutoff>100) return;
  ParamTreeCons * paramtreecons=new ParamTreeCons();
  paramtreecons->setBlocks(outputFile->getBlocks());
  paramtreecons->setNames(outputFile->getNames());

  paramtreecons->setTimeScale(param->getTimeScale());
  if(param) delete(param);
  param=paramtreecons;
  outputFile->startOver();
  while (outputFile->getIt(paramtreecons))
  {
    displayStatus();
    paramtreecons->account();
  }
  paramtreecons->consensus(cutoff);
  param->newDisplayTree(param->getTree(),false);
  if(dest.length()==0) {
	cout<<param->getDisplayTree()->newick()<<endl;// Included for old functionality - do we want it?
  	repaint();
  }else{
	ostream* f_out;
	if(dest.compare("COUT")==0) cout<<param->getDisplayTree()->newick()<<endl;
	else {f_out = new ofstream(qPrintable(dest));
	*f_out<<param->getDisplayTree()->newick()<<endl;
	if(dest.compare("COUT")!=0)delete f_out;}
  }
}

void MainWindowImpl::on_actionExtended_consensus_of_trees_activated(QString dest)
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  ParamTreeCons * paramtreecons=newExtConsTree();
  if(param) delete(param);
  param=paramtreecons;
  if(dest.length()==0) {
	cout<<param->getDisplayTree()->newick()<<endl;// Included for old functionality - do we want it?
  	repaint();
  }else{
	ostream* f_out;
	if(dest.compare("COUT")==0) cout<<param->getDisplayTree()->newick()<<endl;
	else {f_out = new ofstream(qPrintable(dest));
	*f_out<<param->getDisplayTree()->newick()<<endl;
	if(dest.compare("COUT")!=0)delete f_out;}
  }
  param->newDisplayTree(param->getTree(),false);
}

void MainWindowImpl::mouseDoubleClickEvent ( QMouseEvent * event ) 
{
	if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
	statusBar()->showMessage("Capturing node...");
	int id=param->getIdAt(event->x(),event->y(),this);
	displayStatus();
	emit wasDbClicked(id,event);
}

void MainWindowImpl::computeAfterDbClick(int id,QMouseEvent * event)
{
	bool getto=true;
	statusBar()->showMessage("Computing...");
  if(!outputFile->isinitialised()) outputFile->countVectors();
if (event->modifiers()==Qt::CTRL) getto=false;
	PlotImpl * pi=new PlotImpl(this);
	if (event->button() && event->modifiers()==Qt::SHIFT) pi->setValues(outputFile->getRelGenoRec(param,id));
	else pi->setValues(outputFile->getGenoRec(id,getto));
	pi->setBlocks(param->getData()->getBlocks());
	pi->setMode(2);
	if (param->isCons) {
param->unsetDisplayTree();
param->newDisplayTree(param->getTree());
  outputFile->startOver();
  ((ParamCons*)param)->reset();
  while (outputFile->getIt(param))
    {
      displayStatus();
      ((ParamCons*)param)->account(id,getto);
    }
  ((ParamCons*)param)->set();
  repaint();
 }
	pi->show();
	displayStatus();
}

void MainWindowImpl::on_actionColour_density_plot_activated()
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}
  disconnect(this,SIGNAL(wasDbClicked(int,QMouseEvent*)),0,0);
  ColouredImpl * ci=new ColouredImpl(this);
  connect(this,SIGNAL(wasDbClicked(int,QMouseEvent*)),ci,SLOT(addGroup(int)));
  ci->show();
}

ParamTreeCons *  MainWindowImpl::newExtConsTree()
{
	if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return(NULL);}
	ParamTreeCons * paramtreecons=new ParamTreeCons();
	paramtreecons->setBlocks(outputFile->getBlocks());
	paramtreecons->setNames(outputFile->getNames());
	paramtreecons->setTimeScale(param->getTimeScale());
	ParamQt * oldparam=param;
 	param=paramtreecons;
	outputFile->startOver();
	while (outputFile->getIt(paramtreecons))
  	{
	  displayStatus(tr("Computing consensus tree..."));
    	  paramtreecons->account();
  	}
  	paramtreecons->consensusExt();
	paramtreecons->setDisplayTree(paramtreecons->getTree());
	param=oldparam;
	return(paramtreecons);
}

void MainWindowImpl::doColourPlot(QStringList*nodes,QStringList*colors,bool denDep,bool colDep)
{
  ParamConsMult * paramcons=new ParamConsMult(nodes,colors,denDep,colDep);
  paramcons->setBlocks(outputFile->getBlocks());
  paramcons->setNames(outputFile->getNames());
  paramcons->setTimeScale(param->getTimeScale());
  paramcons->setGene(param->getGene());
  if(param->displaySet()) paramcons->newDisplayTree(param->getDisplayTree(),false);
  else paramcons->newDisplayTree(param->getTree());

  if(param) delete(param);
  param=paramcons;
  
  outputFile->startOver();
  while (outputFile->getIt(param))
    {
      displayStatus();
      paramcons->account();
    }
  paramcons->set();
  repaint();
}

void MainWindowImpl::on_actionHeat_map_activated(int correctforprior,QString dest)
{
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}

  if(!param->displaySet()) param->newDisplayTree(param->getTree());
  outputFile->startOver();
  HeatImpl * hi=new HeatImpl(param->getTree()->getN()*2-1);
  while (outputFile->getIt(param))
  {
    displayStatus();
    hi->account(param);
  }
  if(correctforprior>0) hi->compute_correct(correctforprior);
  else hi->compute();
  if(dest.length()==0) {
	hi->show();
  }else {
	ostream* f_out;
	if(dest.compare("COUT")==0) f_out=&cout;
	else f_out = new ofstream(qPrintable(dest));
	hi->print(f_out);
	if(dest.compare("COUT")!=0) delete f_out;
  }
}

void MainWindowImpl::on_actionHeat_map_activated(){ 
on_actionHeat_map_activated(false,QString(""));
}

void MainWindowImpl::on_actionPd_map_activated(){
on_actionPd_map_activated(0,QString(""));
}

void MainWindowImpl::on_actionPheat_map_rel_activated(){
on_actionPheat_map_activated(0,QString(""),true);
}

void MainWindowImpl::on_actionPheat_map_activated(){
on_actionPheat_map_activated(0,QString(""),false);
}

void MainWindowImpl::on_actionPd_map_activated(int correctforprior,QString dest){
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}

  bool ok;
  int sitejump=QInputDialog::getInteger(this,"Site Sample Rate","Enter gap between evaluated sites:",1,1,10000,1,&ok);
  if(!ok) return;
  if(!param->displaySet()) param->newDisplayTree(param->getTree());
  outputFile->startOver();
  PdImpl * pd=new PdImpl(param->getTree()->getN(),sitejump);
 // outputFile->getIt(param);
  while (outputFile->getIt(param))
  {
    displayStatus();
    pd->account(param);
  }
  pd->compute(correctforprior);

  if(dest.length()==0) {
	pd->show();
  }else {
	ostream* f_out;
	if(dest.compare("COUT")==0) f_out=&cout;
	else f_out = new ofstream(qPrintable(dest));
	pd->print(f_out);
	if(dest.compare("COUT")!=0)delete f_out;
  }
}

void MainWindowImpl::on_actionPheat_map_activated(int correctforprior,QString dest,bool reldists){
  if (outputFile==NULL) {QMessageBox::about(0,"Information","Need some data first.");return;}

  if(!param->displaySet()) param->newDisplayTree(param->getTree());
  outputFile->startOver();
  PHeatImpl * pd=new PHeatImpl(param->getTree()->getN(),reldists);
//  outputFile->getIt(param);
  while (outputFile->getIt(param))
  {
    displayStatus();
    pd->account(param);
  }
  pd->compute(correctforprior);

  if(dest.length()==0) {
	pd->show();
  }else {
	ostream* f_out;
	if(dest.compare("COUT")==0) f_out=&cout;
	else f_out = new ofstream(qPrintable(dest));
	pd->print(f_out);
	if(dest.compare("COUT")!=0)delete f_out;
  }
}


void MainWindowImpl::on_actionTest_convergence_of_trees_activated(QStringList files,QString dest)
{
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  if(files.size()==0) files  = QFileDialog::getOpenFileNames(this, tr("Select File(s)"),".","XML files (*.xml);;All files (*)");
  if (files.size()==0) return;
  GelmanRubinImpl * gr=new GelmanRubinImpl(this);
  files.push_front(outputFile->getFileName());

  if(dest.length()==0) {
	gr->computeTree(param,&files);
  	gr->show();
  }else if(dest.compare("COUT")==0) gr->computeTree(param,&files,&cout);
  else {
	ostream* f_out;
	f_out = new ofstream(qPrintable(dest));
	gr->computeTree(param,&files,f_out);
//	*f_out<<qPrintable(qs)<<endl;
	delete f_out;
  }
  loadIteration();
}

void MainWindowImpl::on_actionScore_against_true_tree_activated(QString qstr,QString dest)
{
  if (outputFile==NULL) {QMessageBox::about(0, "Information","Need some data first.");return;}
  if(qstr.length()==0) qstr = QFileDialog::getOpenFileName(this, tr("Select tree file"),".","Newick file (*.nwk *.newick *.tree);;All files (*)");

  if (qstr==NULL) cout<<"qstr was null"<<endl;
  if (qstr==NULL) return;
  QString qs=GelmanRubinImpl::compareTrueTree(param,outputFile->getFileName(),qstr);
  if(dest.length()==0) QMessageBox::about(0,"Result of comparison with true tree",qs);
  else if(dest.compare("COUT")==0) cout<<qPrintable(qs)<<endl;
  else {
	ostream* f_out;
	f_out = new ofstream(qPrintable(dest));
	*f_out<<qPrintable(qs)<<endl;
	delete f_out;
  }
  loadIteration();
}

void MainWindowImpl::displayStatus(QString str)
{
  statusBar()->showMessage(tr("Sample ")+QString::number(outputFile->getCurIt()+1)
  +tr(", Iteration ")+QString::number(param->getNumber())+(", Timescale=")+QString::number(param->getTimeScale())+tr(" ") + str);
}

void MainWindowImpl::on_actionNext_site_activated()
{
  explorerSite++;if (explorerSite>=param->getData()->getL()) explorerSite=param->getData()->getL()-1;
  ((ParamMR*)param)->consensus(explorerCutoff,explorerSite);
  statusBar()->showMessage(tr("Explorer mode, site=")+QString::number(explorerSite));
  repaint();
}

void MainWindowImpl::on_actionPrev_site_activated()
{
  explorerSite--;if (explorerSite<0) explorerSite=0;
  ((ParamMR*)param)->consensus(explorerCutoff,explorerSite);
  statusBar()->showMessage(tr("Explorer mode, site=")+QString::number(explorerSite));
  repaint();
}

void MainWindowImpl::on_actionSet_cutoff_activated()
{
  bool ok;
  int cutoff=QInputDialog::getInteger(this,"Enter cutoff","Enter value of cutoff:",explorerCutoff,0,100,1,&ok);
  if (!ok) return;
  explorerCutoff=cutoff;
  ((ParamMR*)param)->consensus(explorerCutoff,explorerSite);
  repaint();
}

void MainWindowImpl::on_actionJump_to_site_activated()
{
  bool ok;
  int site=QInputDialog::getInteger(this,"Enter site","Enter number of site:",explorerSite,0,param->getData()->getL()-1,1,&ok);
  if (!ok) return;
  jumpToSite(site);
}

void MainWindowImpl::jumpToSite(int site)
{
  explorerSite=site;
  ((ParamMR*)param)->consensus(explorerCutoff,explorerSite);
  statusBar()->showMessage(tr("Explorer mode, site=")+QString::number(explorerSite));
  repaint();
}

void MainWindowImpl::on_actionExport_movie_activated()
{
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save movie file"),".",tr("WMV files (*.wmv);;AVI files (*.avi)"));
  if (qstr==NULL) return;
  bool ok;int step=QInputDialog::getInteger(this,"Enter step","Enter value of step:",1,1,100000,1,&ok);if (!ok) return;
  QImage image(width(),height(),QImage::Format_ARGB32);
  explorerSite=0;
  QString time=QTime::currentTime().toString();
  while (explorerSite<=param->getData()->getL())
  {
  image.fill(0);
  image.invertPixels();//Fill image in white
  statusBar()->showMessage(tr("Explorer mode, site=")+QString::number(explorerSite));
  ((ParamMR*)param)->consensus(explorerCutoff,explorerSite);
  param->display(&image);
  QPainter painter(&image);
  painter.drawText(image.width()*0.05,image.height()*0.05,image.width()*0.1,image.height()*0.1,Qt::AlignCenter,QString::number(explorerSite));
  image.save(QString("/tmp/%1it%2.jpg").arg(time).arg(QString::number(explorerSite/step),7,'0'));
  explorerSite+=step;
  }
  explorerSite=0;
  QProcess ffmpeg;
  if (qstr.endsWith("wmv")) 
  ffmpeg.start("ffmpeg", QString("-y -r 10 -i /tmp/%1it%07d.jpg -vcodec wmv1 %2").arg(time).arg(qstr).split(" "));
  else 
  ffmpeg.start("ffmpeg", QString("-y -r 10 -i /tmp/%1it%07d.jpg %2").arg(time).arg(qstr).split(" "));
  ffmpeg.waitForFinished();
}

void MainWindowImpl::on_actionExport_CSV_activated()
{
  QString qstr = QFileDialog::getSaveFileName(this, tr("Export to CSV"),".","CSV index files (*.csv);;All files (*)");
  if (qstr==NULL) return;
  bool ok;int step=QInputDialog::getInteger(this,"Enter step","Enter value of step:",1,1,100000,1,&ok);if (!ok) return;
  QFile file(qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream out(&file);
  ((ParamMR*)param)->toCSV(&out,explorerCutoff,step);
  file.close();
}

void MainWindowImpl::on_actionExport_to_Artemis_activated()
{
  QString qstr = QFileDialog::getSaveFileName(this, tr("Export to Artemis"),".","GBK files (*.gbk);;All files (*)");
  if (qstr==NULL) return;
  QFile file(qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream out(&file);
  ((ParamMR*)param)->toArtemis(&out,explorerCutoff);
  file.close();
}

void MainWindowImpl::on_actionCompute_prior_corrected_activated()
{
	computeExplorer(true);
}

void MainWindowImpl::on_actionNameFormatm1_activated()
{
	actionNameFormatm1->setChecked(true);
	actionNameFormat0->setChecked(false);
	actionNameFormat1->setChecked(false);
	actionNameFormat2->setChecked(false);
	param->setNameType(-1);
	repaint();
}

void MainWindowImpl::on_actionNameFormat0_activated()
{
	actionNameFormatm1->setChecked(false);
	actionNameFormat0->setChecked(true);
	actionNameFormat1->setChecked(false);
	actionNameFormat2->setChecked(false);
	param->setNameType(0);
	repaint();
}

void MainWindowImpl::on_actionNameFormat1_activated()
{
	actionNameFormatm1->setChecked(false);
	actionNameFormat0->setChecked(false);
	actionNameFormat1->setChecked(true);
	actionNameFormat2->setChecked(false);
	param->setNameType(1);
	repaint();
}

void MainWindowImpl::on_actionNameFormat2_activated()
{
	actionNameFormatm1->setChecked(false);
	actionNameFormat0->setChecked(false);
	actionNameFormat1->setChecked(false);
	actionNameFormat2->setChecked(true);
	param->setNameType(2);
	repaint();
}
