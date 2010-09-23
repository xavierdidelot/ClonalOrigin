#include "gelmanrubinimpl.h"
//
GelmanRubinImpl::GelmanRubinImpl( QWidget * parent, Qt::WFlags f)
    : QDialog(parent, f) {
  setupUi(this);
  param=NULL;
  outputFile=NULL;
  out=NULL;
  sep=",";
  doGR=true;
  groupBox->setEnabled(true);
//    groupBox->setMinimumSize(QSize(0, 0));
    groupBox->setHidden(true);
 //connect(pushGo,SIGNAL(wasDbClicked(int,QMouseEvent*)),this,SLOT(computeAfterClick(int,QMouseEv
	//QObject::connect(pushGo, SIGNAL(clicked()), auxSignals, SLOT(buttonClicked()));
}
//

void GelmanRubinImpl::compute(ParamQt*param,OutputFile*outputFile,QStringList*others,ostream* out,bool getparams,bool getnumedges,bool getpairwisedists) {
  vector< vector< vector<double> > > v;//Store:file x iteration x property
  int m=1+others->size();//Number of files

  v.push_back(vector< vector<double> >());
  outputFile->startOver();
  while (outputFile->getIt(param)) {
	v[0].push_back(extractInfo(param,false,getparams,getnumedges,getpairwisedists));
  }
  for (int i=0;i<others->size();i++) {
      OutputFile * out2=new OutputFile(QStringList(others->at(i)),false);
      v.push_back(vector< vector<double> >());
      while (out2->getIt(param)) {
	v[i+1].push_back(extractInfo(param,false,getparams,getnumedges,getpairwisedists));
      }
    }

  //Table headers
  table->setRowCount(1+v[0][0].size());
  table->setColumnCount(m+2);
  table->setColumnWidth(0,table->columnWidth(0)*2);
  table->setItem(0,0,new QTableWidgetItem(QString("Parameter")));
  for (int i=1;i<=m;i++) table->setItem(0,i,new QTableWidgetItem(QString("File ")+QString::number(i)));
  table->setItem(0,1+m,new QTableWidgetItem(QString("Gelman-Rubin")));

  nameTable(param,getparams,getnumedges,getpairwisedists);
  //Table
  for (unsigned int whichP=0;whichP<v[0][0].size();whichP++) {
      vector< vector<double> > data;
      for (int i=0;i<m;i++) {
          data.push_back(vector<double>());
          for (unsigned int j=0;j<v[i].size();j++) {
              data[i].push_back(v[i][j][whichP]);
            }
        }
      double R=test(&data);
      for (int i=1;i<=m;i++) {
          double mean=0.0;
          for (unsigned int j=0;j<data[i-1].size();j++) mean+=data[i-1][j]/data[i-1].size();
          table->setItem(whichP+1,i,new QTableWidgetItem(QString::number(mean)));
        }
      table->setItem(whichP+1,1+m,new QTableWidgetItem(QString::number(R)));
    }
  if(out!=NULL){
	for(int i=0;i<table->rowCount();i++){
	  for(int j=0;j<table->columnCount()-1;j++){
		*out<<table->item(i,j)->text().toStdString()<<sep;
	  }
	*out<<table->item(i,table->columnCount()-1)->text().toStdString()<<endl;
	}
  }
}

void GelmanRubinImpl::computeTree(ParamQt*param,QStringList*files,ostream* out) {
  int m=files->size();//Number of files
  int iso=param->getTree()->getN();//Number of isolates
  QHash<QString, vector<double> > hash;
  for (int i=0;i<files->size();i++) {
      OutputFile * out2=new OutputFile(QStringList(files->at(i)),false);
      int its=0;
      while (out2->getIt(param)) {
          its++;
          for (int ii=iso;ii<iso*2-2;ii++) {
              Node * node=param->getTree()->getNode(ii);
              QString key(iso,'0');
              ParamTreeCons::makeKey(node,&key);
              if (!hash.contains(key)) hash[key]=vector<double>(m,0.0);
              hash[key][i]++;
            }
        };
      QStringList keys=hash.uniqueKeys();
      for (int ii=0;ii<keys.size();ii++) hash[keys[ii]][i]/=its;
    }

  vector< vector< vector<double> > > v;//Store:file x iteration x property
  v.push_back(vector< vector<double> >());

  QStringList keys=hash.uniqueKeys();

  //Average standard deviation of split frequencies
  double asdsf=0.0;
  int cmp=0;
for (int i=0;i<m;i++) for (int j=i+1;j<m;j++) for (int k=0;k<keys.size();k++) {asdsf+=gsl_pow_2(hash[keys[k]][i]-hash[keys[k]][j]);cmp++;}
  asdsf=sqrt(asdsf/cmp);
  label->setText(QString("Average standard deviation of split frequencies=")+QString::number(asdsf));

  //Table headers
  table->setRowCount(keys.size()+1);
  table->setColumnCount(m+1);
  table->setColumnWidth(0,table->columnWidth(0)*2);
  if(out!=NULL) *out<<asdsf<<",";
  for (int i=1;i<=m;i++) {
	table->setItem(0,i,new QTableWidgetItem(QString("File ")+QString::number(i)));
	if(out!=NULL) *out<<"File "<<i;
	if(out!=NULL && i<m) *out<<",";
  }
  if(out!=NULL) *out<<endl;
  for (int i=0;i<keys.size();i++)
    table->setItem(i+1,0,new QTableWidgetItem(keys[i]));

  //Table
  for (int i=0;i<keys.size();i++){
    if(out!=NULL) *out<<qPrintable(keys[i])<<",";
    for (int j=0;j<m;j++) {
	table->setItem(i+1,j+1,new QTableWidgetItem(QString::number(hash[keys[i]][j])));
	  if(out!=NULL) *out<<hash[keys[i]][j];
	  if(out!=NULL && j<m-1) *out<<",";
    }
    if(out!=NULL) *out<<endl;
  }
}

void GelmanRubinImpl::nameTable(ParamQt*p,bool getparams,bool getnumedges,bool getpairwisedists) {
  int rownumber=1;
  int iso=p->getTree()->getN();//Number of isolates
  if(getparams){
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("theta")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("rho")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("delta")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("Num Edges")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("Likelihood")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("TMRCA")));
  table->setItem(rownumber++,0,new QTableWidgetItem(QString("SumBraLen")));
  }
  if(getnumedges){
    for (int i=0;i<iso+iso-2;i++)
      table->setItem(rownumber++,0,new QTableWidgetItem(QString("Edges on branch ")+QString::number(i)));
  }
  if(getpairwisedists){
	for (int i=0;i<iso;i++)for (int j=i+1;j<iso;j++)
      table->setItem(rownumber++,0,new QTableWidgetItem(QString("PairwiseDistance")+QString::number(i)+QString("to")+QString::number(j)));
  }
}


vector<double> GelmanRubinImpl::extractInfo(ParamQt*p, bool csv,bool getparams,bool getnumedges,bool getpairwisedists) {
  vector<double> v;
  int ccount=0;
  if (csv) ccount=p->countConv();
  int iso=p->getTree()->getN();
  if(getparams) { for (int whichP=0;whichP<7+ccount;whichP++) {
      if (whichP==0) {v.push_back(p->getTheta());continue;}
      if (whichP==1) {v.push_back(p->getRho  ());continue;}
      if (whichP==2) {v.push_back(p->getDelta());continue;}
      if (whichP==3) {v.push_back(p->getTree()->numRecEdge());continue;}
      if (whichP==4) {v.push_back(p->getLL());continue;}
      if (whichP==5) {v.push_back(p->getTree()->getNode(iso+iso-2)->getAge());continue;}
      if (whichP==6) {v.push_back(p->getTree()->getTTotal());continue;}
      if (whichP>6 && whichP<7+ccount) {v.push_back(p->getConvData(whichP-7));continue;}
  } }
  if(getnumedges) {for (int whichP=0;whichP<(iso+iso-2);whichP++) {
       v.push_back(p->getTree()->numRecEdgeOnBranch(whichP));
  }}
  if(getpairwisedists) {
	vector<double> vt=p->pairwiseDistanceList();
	for(unsigned int c1=0;c1<vt.size();c1++) v.push_back(vt[c1]);
  }
  if(!getparams &&!getnumedges && !getpairwisedists) {
	vector<double> vt;
        p->consistentAgeList(&vt);
	for(unsigned int c1=0;c1<vt.size();c1++) v.push_back(vt[c1]);
  }
  return v;
}


double GelmanRubinImpl::test(vector< vector<double> >*data) {
  vector<double> means;
  double mean;
  unsigned int m=data->size();
  unsigned int n=data->at(0).size();
 for (unsigned int i=0;i<m;i++) if(data->at(i).size()<n) n=data->at(i).size();

  for (unsigned int i=0;i<m;i++) {
      mean=0.0;
      for (unsigned int j=0;j<n;j++) {
	mean+=data->at(i).at(j);	
      }
      mean/=n;
      means.push_back(mean);
    }
  mean=0.0;
  for (unsigned int j=0;j<m;j++) mean+=means[j];
  mean/=m;
  double B=0.0;
  for (unsigned int j=0;j<m;j++) B+=pow((means[j]-mean),2.0);
  B=B*n/(m-1.0);
  double W=0.0;
  for (unsigned int j=0;j<m;j++) for (unsigned int t=0;t<n;t++) W+=pow((data->at(j).at(t)-means[j]),2.0);
  W=W/m/(n-1.0);
  double V=(n-1.0)/n*W+B/n;
  return V/W;
}

void GelmanRubinImpl::on_exportButton_clicked() {
  QString qstr = QFileDialog::getSaveFileName(this, tr("Save output file"),".","CSV files (*.csv);;All files (*)");
  if (qstr==NULL) return;
  QFile file(qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream ts( &file );
  for (int i=0;i<table->rowCount();i++)
    for (int j=0;j<table->columnCount();j++) {
        if (table->item(i,j)!=NULL) ts<<table->item(i,j)->text();
        if (j<table->columnCount()-1) ts<<","; else ts<<endl;
      }
  file.close();
}

void GelmanRubinImpl::outputTracer(ParamQt*param,OutputFile*outputfile,QString*qstr,bool csv,bool getparams,bool getnumedges,bool getpairwisedists) {
  //Extract the information
  vector< vector<double> > v;//Store:iteration x property
  int iso=param->getTree()->getN();//Number of isolates
  outputfile->startOver();
  while (outputfile->getIt(param)) {
	v.push_back(extractInfo(param,csv,getparams,getnumedges,getpairwisedists));
  }

  //Write header
  int ccount=0;
  if (csv) ccount=param->countConv();
  QFile file(*qstr);
  if ( !file.open(QIODevice::WriteOnly)) return;
  QTextStream ts( &file );
  QString sep;
  if (csv) sep=",";else sep="\t";
  if (csv && getparams) {
      ts<<"\"iter\",\"theta\",\"rho\",\"delta\",\"numrecedge\",\"likelihood\",\"TMRCA\",\"SumBraLen\",";
      for (int i=0;i<ccount;i++) {
	ts<< "\""<< param->getConvName(i).c_str()<< "\"";
	if (i<ccount-1 || getpairwisedists || getparams) ts<<sep;
      }
    } else if(getparams){
	ts<<"iter\ttheta\trho\tdelta\tnumrecedge\tlikelihood\ttmrca\tsumbralen";
	if(getnumedges||getpairwisedists)ts<<"\t";
    }
  if(getnumedges) {for (int i=0;i<iso+iso-2;i++) {
      if (csv) ts<<"\"branch"<<i<<"\""; else ts<<"branch"<<i;
      if (i<iso+iso-3 || getpairwisedists) ts<<sep;
  }}
  if(getpairwisedists) {
	for (int i=0;i<iso;i++) {
	  for (int j=i+1;j<iso;j++) {
		if (csv) ts<<"\"PairwiseDistance"<<i<<"to"<<j<<"\""; else ts<<"PairwiseDistance"<<i<<"to"<<j;
      		if (i<iso+iso-3 || getpairwisedists) ts<<sep;
	  }
	}
  }
  if(!getpairwisedists&& !getnumedges && !getparams) {
	ts<<"\"iter\""<<sep;
	for (int i=0;i<iso-1;i++) {
	  ts<<"\"blen"<<i<<"\"";
	  if (i<iso-2) ts<<sep;
	}
  }
  ts<<endl;

  //Write the information
  int state=0;
  for (unsigned int i=0;i<v.size();i++) {
      state++;
      ts<<state<<sep;
      for (unsigned int j=0;j<v[i].size();j++) {
          ts<<v[i][j];
          if (j<v[i].size()-1) ts<<sep; else ts<<endl;
        }
    }
  file.close();
}

QString GelmanRubinImpl::compareTrueTree(ParamQt*param,QString outputfile,QString newick) {
  int iso=param->getTree()->getN();//Number of isolates
  QHash<QString, vector<double> > hash;
  //Analyze partitions in the output file
  OutputFile * out2=new OutputFile(QStringList(outputfile),false);
  int its=0;
  while (out2->getIt(param)) {
      its++;
      for (int ii=iso;ii<iso*2-2;ii++) {
          Node * node=param->getTree()->getNode(ii);
          QString key(iso,'0');
          ParamTreeCons::makeKey(node,&key);
          if (!hash.contains(key)) hash[key]=vector<double>(2,0.0);
          hash[key][0]++;
        }
    };
  QStringList keys=hash.uniqueKeys();
  for (int ii=0;ii<keys.size();ii++) hash[keys[ii]][0]/=its;

  //Analyze partitions in true tree
  Tree * truth=new Tree(newick.toStdString());
  for (int ii=iso;ii<iso*2-2;ii++) {
      Node * node=truth->getNode(ii);
      QString key(iso,'0');
      ParamTreeCons::makeKey(node,&key);
      if (!hash.contains(key)) hash[key]=vector<double>(2,0.0);
      hash[key][1]=1.0;
    }

  //Calculate stats

  keys=hash.uniqueKeys();
  double score=0.0;
  double efficiency=0.0;double eff2=0.0;
  double accuracy=0.0;double acc2=0.0;
  for (int ii=0;ii<keys.size();ii++) {
  if (hash[keys[ii]][1]==0.0) score-=hash[keys[ii]][0];else score+=hash[keys[ii]][0];
  if (hash[keys[ii]][1]==1.0) {eff2++;if (hash[keys[ii]][0]>=0.5) efficiency++;};
  if (hash[keys[ii]][0]>=0.5) {acc2++;if (hash[keys[ii]][1]==1.0) accuracy  ++;};
  }
  return QString("Score=%1; Efficiency=%2; Accuracy=%3").arg(score).arg(efficiency/eff2).arg(accuracy/acc2);
}


void GelmanRubinImpl::on_buttonBox_accepted() {
  if(!doGR) {
	bool usecsv=true;
	if(sep.compare(",")!=0)usecsv=false;
    if(param==NULL) throw("Export not initalised!");
    if(radioType1->isChecked()) outputTracer(param,outputFile,&(others[0]),usecsv,true,false,false);
    else if(radioType2->isChecked()) outputTracer(param,outputFile,&(others[0]),usecsv,true,true,false);
    else if(radioType3->isChecked()) outputTracer(param,outputFile,&(others[0]),usecsv,true,false,true);
    else if(radioType4->isChecked()) outputTracer(param,outputFile,&(others[0]),usecsv,true,true,true);
    else if(radioType5->isChecked()) outputTracer(param,outputFile,&(others[0]),usecsv,false,false,false);
  }
  close();
}

void GelmanRubinImpl::on_buttonBox_rejected() {
  close();
}

void GelmanRubinImpl::on_pushGo_clicked()
{
    if(param==NULL) throw("Gelman Rubin not initalised!");
    if(radioType1->isChecked()) compute(param,outputFile,&others,NULL,true,false,false);
    else if(radioType2->isChecked()) compute(param,outputFile,&others,NULL,true,true,false);
    else if(radioType3->isChecked()) compute(param,outputFile,&others,NULL,true,false,true);
    else compute(param,outputFile,&others,NULL,true,true,true);
    groupBox->setEnabled(false);

}
