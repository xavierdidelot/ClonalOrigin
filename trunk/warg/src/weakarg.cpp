#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "rectree.h"
#include "param.h"
#include "move.h"
#include "data.h"
#include "rng.h"
#include <cstring>
#include <fstream>
#include "mpiutils.h"
#include "weakarg.h"
#include "wargxml.h"

using namespace std;
namespace weakarg
{
bool initializeTree(Data* &datap, RecTree* &rectree,vector<string> inputfiles,string datafile);///< initialises the rectree based on the inputfiles specified.
RecTree * makeGreedyTree(Data * data,WargXml * infile,vector< vector<double> >  * sumdetails,int *count,vector<double> * pars,vector<double> *sumdists);///< makes a greeedy "best tree" from a previous warg run
vector<double> readInputFiles(Data* &data, RecTree* &rectree,vector<double> &sumdists,vector<int> &keepregions,vector<int> &previousL,vector<string> inputfiles,string datafile,int greedystage);///< Reads the input files and processes them according to the stage of the input procedure.

ProgramOptions& opt() {
	static ProgramOptions po;	// define a single instance of ProgramOptions per process.
	return po;
}


} // end namespace weakarg

string getVersion(){
	string ret;
#ifdef PACKAGE_STRING
	ret.append(PACKAGE_STRING);
#else
	ret.append("warg");
#endif
	ret.append(" build date "); 
	ret.append(__DATE__);
	ret.append(" at ");
	ret.append(__TIME__);
	return(ret);
}

void printVersion(){
	cout<<getVersion()<<endl;
}

using namespace weakarg;
// main function goes outside the weakarg namespace so that
// the linker can find it

static const char * help=
    "\
    Usage: weakarg [OPTIONS] treefile datafile outputfile\n\
    \n\
    Options:\n\
    -w NUM      	Sets the number of pre burn-in iterations (default is 100000)\n\
    -x NUM      	Sets the number of burn-in iterations (default is 100000)\n\
    -y NUM      	Sets the number of iterations after burn-in (default is 100000)\n\
    -z NUM      	Sets the number of iterations between samples (default is 100)\n\
    -T NUM      	Sets the value of theta. Use sNUM instead of NUM for per-site\n\
    -R NUM      	Sets the value of rho. Use sNUM instead of NUM for per-site\n\
    -D NUM      	Sets the value of delta\n\
    -s NUM      	Use given seed to initiate random number generator\n\
    -S NUM,SEED 	Run on a subset of NUM regions determined by seed SEED\n\
       NUM/NUM/../NUM 	Run on a specified region(s) given by each NUM.\n\
    -r NUM		Perform r tempered steps between topological updates (default:0)\n\
    -t NUM		Tempered at \"temperature\" t for topological updates (default:1.0)\n\
    -U			Start from UPGMA tree, rather than the default random tree.\n\
    -G NUM		Greedily compute the \"best fit\" tree, given the recombination\n\
			observed on the current tree.  If NUM is negative and a previous\n\
			run is provided, the tree is calculated from all observed values.\n\
			If NUM is positive, a \"greedy move\" is performed with weight\n\
			NUM (see -a).  Note that this is NOT an MCMC move and causes bias.\n\
    -a NUM,...,NUM	Set the ELEVEN (real valued) move weightings to the given vector,\n\
    with weightings separated by commas (NOT SPACES).  \n\
    The weightings need not sum to 1, but must be in the following order:\n\
    	MoveRho   (ignored if not needed)\n\
    	MoveDelta (ignored if not needed)\n\
    	MoveTheta (ignored if not needed)\n\
    	MoveRemEdge\n\
    	MoveAddEdge\n\
    	MoveSiteChange\n\
    	MoveTimeChange\n\
    	MoveEdgeChange\n\
    	MoveAgeClonal\n\
    	MoveScaleTree\n\
    	MoveRegraftClonal\n\
    -i NUM,...,NUM	Set the SIX parameters for creating random Recombination Trees\n\
			under the inference model.  The parameters are:\n\
    	N	(integer)	The number of sequences in the sample (default 10)\n\
    	n_B	(integer)	The number of block boundaries in the sample (default 8)\n\
    	l_B	(integer)	The length of each block: L=n_B * l_B (default 500)\n\
    	delta	(real)		The average length of imports (default 500.0)\n\
    	theta	(real)		The mutation rate NOT per site (default 100.0)\n\
    	rho	(real)		The recombination rate NOT per site (default 50.0)\n\
    -f			Forbid topology changes, (allowing updates of coalescence times).\n\
    -v          	Verbose mode\n\
    -h          	This help message\n\
    -V          	Print Version info\n\
    ";

int main(int argc, char *argv[])
{
    string comment="Command line: ";
    for(int c1=0;c1< argc;c1++) {comment.append(argv[c1]);comment.append(" ");}
    comment.append("\nVersion: ");
    comment.append(getVersion());
    vector<string> inputfiles;
    initmpi(argc,argv);
	makerng(true);
    optind=0;
    bool upgma=false;
    int c;
    char * pch;
    double simparrho=50.0;
    double simpartheta=100.0;
    double simpardelta=500.0;
    int simparN=10;
    int simparnumblocks=8;
    int simparblocksize=500;
    std::stringstream ss;
    unsigned long seed=0;
    bool readparams=false;
    bool setregions=false;
    while ((c = getopt (argc, argv, "w:x:y:z:s:va:T:R:D:L:C:r:t:i:S:G:fUhV")) != -1)
        switch (c)
        {
        case('w'):if(atoi(optarg)>=0)opt().preburnin=atoi(optarg);break;
        case('x'):if(atoi(optarg)>=0)opt().burnin=atoi(optarg);break;
        case('y'):if(atoi(optarg)>=0)opt().additional=atoi(optarg);break;
        case('z'):if(atoi(optarg)> 0)opt().thinin=atoi(optarg);break;
        case('T'):opt().theta=atof(optarg);if (optarg[0]=='s') {opt().theta=atof(optarg+1);opt().thetaPerSite=true;};break;
        case('R'):opt().rho=atof(optarg);if (optarg[0]=='s') {opt().rho=atof(optarg+1);opt().rhoPerSite=true;};break;
        case('D'):opt().delta=atof(optarg);break;
        case('s'):seed=strtoul(optarg,NULL,10);break;
        case('v'):opt().verbose=true;break;
	case('f'):opt().allowclonal=false;break;
        case('U'):upgma=true;break;
        case('a'):pch = strtok (optarg,",");
            for(int i=0;i<NUMMOVES;i++) {
		if(pch==NULL) {cout<<"Wrong -a string."<<endl<<help<<endl;return 1;}
                opt().movep[i]=fabs(atof(pch));
                pch = strtok (NULL, ",");
            };break;
        case('i'):pch = strtok (optarg,",");
            for(int i=0;i<6;i++) {
		if(pch==NULL) {cout<<"Wrong -i string."<<endl<<help<<endl;return 1;}
		switch(i){
                	case(0):simparN=atoi(pch);break;
                	case(1):simparnumblocks=atoi(pch);break;
                	case(2):simparblocksize=atoi(pch);break;
                	case(3):simpardelta=fabs(atof(pch));break;
                	case(4):simpartheta=fabs(atof(pch));break;
                	case(5):simparrho=fabs(atof(pch));break;
			case '?':cout<<"Wrong -i string."<<endl<<help<<endl;return 1;
		}
                pch = strtok (NULL, ",");
            };break;
	case('r'):opt().temperreps=atoi(optarg);break;
	case('t'):opt().temperT=atof(optarg);break;
        case('L'):opt().logfile=optarg;break;
        case('C'):opt().csvoutfile=optarg;break;
        case('V'):printVersion();  return 0;
	case('S'):pch= strrchr (optarg,',');
	if(pch!=NULL){pch = strtok (optarg,",");opt().subset.push_back(atoi(pch));pch = strtok (NULL,",");opt().subsetSeed=atoi(pch);
	}else{
	    pch = strtok (optarg,"/");
	    while (pch != NULL) {
		opt().subset.push_back(atoi(pch));
    		pch = strtok (NULL, "/");
  	    }
	}
	setregions=true;
	break;
	case('G'):opt().greedyWeight=atof(optarg);break;
	case('h'):cout<<help<<endl;return 0;
        case '?':cout<<"Wrong arguments: did not recognise "<<c<<" "<<optarg<<endl<<help<<endl;return 1;
        default:
            abort ();
        }
    seed=seedrng(seed);// <0 means use /dev/random or clock.
    comment.append("\nSeed: ");
    ss<<seed;
    comment.append(ss.str());
    if (argc-optind==3 || (opt().greedyWeight<0 && argc-optind>3)){
	while(argc-optind>2) inputfiles.push_back(string(argv[optind++]));
    }
if (argc-optind!=1 && argc-optind!=2) {cout<<"Wrong number of arguments."<<endl<<help<<endl;return 1;}

    Param p;
    RecTree*rectree=NULL;
    Data*data=NULL;
    if (argc-optind==1) {//Run on simulated tree and data
        dlog(1)<<"Simulating rectree..."<<endl;
        vector<int> blocks;
        for (int i=0;i<simparnumblocks;i++) blocks.push_back(i*simparblocksize);
        rectree=new RecTree(simparN,simparrho,simpardelta,blocks);
        dlog(1)<<"Initiating parameter"<<endl;
        p=Param(rectree,NULL);
        dlog(1)<<"Simulating data..."<<endl;
        p.setTheta(simpartheta);
        p.simulateData(blocks);
        p.setTheta(-1.0);
        data=p.getData();
        ofstream dat;
        dat.open("simulatedData.xmfa");
        data->output(&dat);
        dat.close();
        ofstream tru;
        tru.open("truth.xml");
	p.setRho(simparrho);
	p.setTheta(simpartheta);
 	p.exportXMLbegin(tru,comment);
        //p.exportXMLbegin(tru);
        p.exportXMLiter(tru);
        p.exportXMLend(tru);
        tru.close();
	// alternative initialisation options
        //while (p.getRecTree()->numRecEdge()>0) p.getRecTree()->remRecEdge(i);// blank tree
	//for(int i=0;i<p.getRecTree()->numRecEdge();i++) {// remove all but a specific edge	if(p.getRecTree()->getRecEdge(i)->getEdgeTo()!=28) {p.getRecTree()->remRecEdge(i);i--;} }
    } else if(argc-optind==2 && inputfiles.size()==0){//Load data from files
        string datafile=string(argv[optind++]);
        dlog(1)<<"Loading data..."<<endl;
        try{data=new Data(datafile);
	}catch(const char *){exit(0);}
        if (upgma)
        {
        rectree=new RecTree(data,0.0,500.0,*(data->getBlocks()));
	data->subset(opt().subset,opt().subsetSeed);
        }else {
	data->subset(opt().subset,opt().subsetSeed);
        dlog(1)<<"Creating random tree..."<<endl;
        rectree=new RecTree(data->getN(),0.0,500.0,*(data->getBlocks()));
        }
        dlog(1)<<"Initiating parameter..."<<endl;
        p=Param(rectree,data);
	p.setRho(0);
    }else{//Load tree and data from files
        string datafile=string(argv[optind++]);
try{
	readparams=initializeTree(data,rectree,inputfiles,datafile);// initialises data and rectree! (passed by reference)
}catch(char *x){cout<<x<<endl;}
	if(data==NULL) {cerr<<"Error: No Data initialised.  Was there a problem with the input file?"<<endl; exit(0);}
	if(rectree==NULL) {cerr<<"Error: No Rectree initialised.  Was there a problem with the input file?"<<endl; exit(0);}
	dlog(1)<<"Initiating parameter..."<<endl;
	p=Param(rectree,data);
	p.setRho(0);
	if(readparams) {
		p.readProgramOptions();
		WargXml infile(inputfiles[0]);
		p.readParamsFromFile(&infile);
	}
    }
    opt().outfile = argv[optind++];

    if(opt().preburnin>0 && (opt().movep[8]>0 || opt().movep[10]>0)) {
	cout<<"Starting Pre-burnin Metropolis-Hastings algorithm.."<<endl;
	double rho=opt().rho;
	opt().rho=0;
	long int burnin= opt().burnin, additional=opt().additional,temperreps=opt().temperreps;
	opt().burnin=opt().preburnin;
	opt().additional=0;
	p.readProgramOptions();
    	p.metropolis(comment);
	opt().burnin=burnin;
	opt().additional=additional;
	opt().rho=rho;
	opt().temperreps=temperreps;
	p.readProgramOptions();
    }else if(!readparams) p.readProgramOptions();
    cout<<"Starting Metropolis-Hastings algorithm............."<<endl;
    p.metropolis(comment);

    dlog(1)<<"Cleaning up..."<<endl;
    if(p.getRecTree()) delete(p.getRecTree());
    if(data) delete(data);
    gsl_rng_free(rng);

    endmpi();
    return 0;
}


namespace weakarg
{

RecTree * makeGreedyTree(Data * data,WargXml * infile,vector< vector<double> >  * sumdetails,int *count,vector<double> * pars,vector<double> *sumdists)
{
vector< vector<double> > tmpmut;
	RecTree * rectree=NULL;
	Param *p=NULL;
	infile->restart();
	std::streampos sp=infile->tellg(),lastsp=sp;
	for(unsigned int i=0;i<sumdetails->size();i++) for(unsigned int j= 0;j<sumdetails->at(i).size();j++)sumdetails->at(i)[j]*=(double)(*count);
	for(unsigned int i=0;i<pars->size();i++) pars->at(i)*=(double)(*count);
	while(!infile->eof() && sp>=0) {
		sp=infile->gotoLineContaining("<Iteration>",false);
		if(infile->eof() || sp<0) break;
		infile->seekg(sp);
		if(rectree!=NULL) delete(rectree);
try{
		rectree=new RecTree(data->getL(),infile);
		lastsp=sp;
		if(p!=NULL) delete(p);
        	p= new Param(rectree,data);
		p->setRho(0);
		p->readProgramOptions();
		p->readParamsFromFile(infile,sp);
		vector< vector<double> > * mutpairwise=p->greedyPairwiseDetails();
		pars->at(0) += p->empiricalRho();
		pars->at(1) += p->empiricalDelta();
		pars->at(2) += p->empiricalTheta(mutpairwise);

		if(sumdetails->size()==0) {
			for(unsigned int i=0;i<mutpairwise->size();i++) sumdetails->push_back(mutpairwise->at(i));
			(*count)=1;
		}else {
			for(unsigned int i=0;i<sumdetails->size();i++) for(unsigned int j=0;j<sumdetails->at(i).size();j++) sumdetails->at(i)[j]+=mutpairwise->at(i)[j];
			(*count)++;
		}
		sp=infile->gotoLineContaining("</Iteration>",false);
		infile->seekg(sp);
}catch(char * x){cerr<<"Error making greedy tree: "<<x<<endl;exit(0);}
	}
	for(unsigned int i=0;i<sumdetails->size();i++) for(unsigned int j= 0;j<sumdetails->at(i).size();j++)sumdetails->at(i)[j]/=(double)(*count);
	for(unsigned int i=0;i<pars->size();i++) pars->at(i)/=(double)(*count);

	*sumdists = vector<double>(p->greedyCalcDists(sumdetails->at(1),sumdetails->at(0)));
	p->greedyApply(*sumdists);
	delete(p);
	infile->clear();
	infile->seekg(lastsp);
	return(rectree);
}


vector<double> readInputFiles(Data* &data, RecTree* &rectree,vector<double> &sumdists,vector<int> &keepregions,vector<int> &previousL,vector<string> inputfiles,string datafile,int greedystage)
{
	bool setregions=false;
	if(opt().subset.size()>0 || opt().subsetSeed !=-1) setregions=true;
	vector <vector<double> >sumdetails;
	vector<double>pars(3,0.0); // parameters
	int counts=0;// counts for the parameters

	dlog(1)<<"Loading data: "<<datafile<<endl;
	for(unsigned int c1=0;c1<inputfiles.size();c1++) {
		if(data!=NULL){ delete(data);}
		data=new Data(datafile); // we have to keep reloading the data
		dlog(1)<<"Loading tree "<<c1<<"... "<<inputfiles[c1]<<endl;
		string treefile=inputfiles[c1];
		WargXml infile(treefile);
		if(infile.isempty()) {cerr<<"Warning: file "<<treefile<<" is empty. Skipping."<<endl;continue;}
		if(infile.gotoLineContaining("<Iteration>",true)<0) {// is a newick file
			if(inputfiles.size()>1) {cerr<<"Warning: multiple newick files given.  Only the final one will be used"<<endl;}
			data->subset(opt().subset,opt().subsetSeed);// apply the subset as provided on the command line
			rectree=new RecTree(data->getL(),treefile);
		}else{// is an xml output file
			if(opt().subset.size()>0) data->subset(opt().subset,opt().subsetSeed);
			else data->readRegionsFromFile(&infile);
			if(greedystage!=2){ // second pass
				previousL.push_back(previousL.back()+data->getL());
				for(unsigned int c2=0;c2<data->getRegions()->size();c2++) 	keepregions.push_back(data->getRegions()->at(c2));
			}
			if(greedystage==0) {// not greedy
				if(rectree!=NULL) delete(rectree);
				rectree=new RecTree(data->getL(),&infile);
			}else if(greedystage==1){// get the dists for a greedy tree
				if(rectree!=NULL) delete(rectree);
				rectree=makeGreedyTree(data,&infile,&sumdetails,&counts,&pars,&sumdists);
			}else if(greedystage==2){// construct a final iteration from all input files
				if(rectree!=NULL && c1==0) delete(rectree);
				if(c1==0) rectree=new RecTree(previousL.back(),&infile,false);
				rectree->addEdgesFromFile(&infile,previousL[c1]);
			}
		}
	}
	return(pars);
}

bool initializeTree(Data* &data, RecTree* &rectree,vector<string> inputfiles,string datafile)
{
	vector<double>pars(3,0.0);
	vector<double> sumdists;
	vector<int> keepregions;
	vector<int> previousL(1,0);// List of partial L's; starts with just 0
	bool  readparams=false;

	if(opt().greedyWeight<0) {// create a greedy tree from the input
	  pars=readInputFiles(data,rectree,sumdists,keepregions,previousL,inputfiles,datafile,1);
	  readInputFiles(data,rectree,sumdists,keepregions,previousL,inputfiles,datafile,2);
	}else{// just read in the input and keep the specified regions
	if(inputfiles.size()>1) cerr<<"Warning: multiple input files specified but this is only purposeful with the -G -1 option. Ignoring all but the final one."<<endl;
	  pars=readInputFiles(data,rectree,sumdists,keepregions,previousL,inputfiles,datafile,0);
	}
	for(unsigned int i=0;i<pars.size();i++) if(pars[i]!=0) readparams=true;

	if(data!=NULL){ delete(data);}
	data=new Data(datafile); // we have to keep reloading the data
	data->subset(keepregions,-1);// apply the subset of all data we've seen
        Param * p= new Param(rectree,data);
	if(readparams) p->setTheta(pars[2]);p->setRho(pars[0]);p->setDelta(pars[2]);
	if(opt().greedyWeight<0) p->greedyApply(sumdists);
	delete(p);
	return(readparams);
}


}
