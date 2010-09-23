#include "rectree.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "slotallocator.h"
#include "mpiutils.h"
#include "param.h"


//#define DEBUG

using namespace std;
namespace weakarg
{

inline double fround(double n, double d)
{
    return floor(n * pow(10., d) + .5) / pow(10., d);
}

RecTree::RecTree(unsigned int numsites,string newick,bool isFilename,bool forceages):Tree(newick,isFilename,forceages)
{
    tabNode=vector<int>(n*2-1);
    tabRec=vector<int>(MAXNODES);
    tabSons=vector<vector<int> >(MAXNODES,vector<int>(2,0));
    tabSonsDist=vector<vector<double> >(MAXNODES,vector<double>(2,0.0));
    tabFather=vector<int>(MAXNODES,-1);
    age=vector<double>(MAXNODES,0.0);
    tabEdge=vector<int>(MAXNODES);
    L = numsites;
    sameLTasPrev=vector<bool>(L,true);
    edge=vector<RecEdge*>(0);
}

RecTree::RecTree(unsigned int numsites,WargXml *infile,bool addedges,bool forceages,bool isfinal,std::streampos itstart):Tree(infile->getTree(itstart,isfinal),false,forceages)
{
    tabNode=vector<int>(n*2-1);
    tabRec=vector<int>(MAXNODES);
    tabSons=vector<vector<int> >(MAXNODES,vector<int>(2,0));
    tabSonsDist=vector<vector<double> >(MAXNODES,vector<double>(2,0.0));
    tabFather=vector<int>(MAXNODES,-1);
    age=vector<double>(MAXNODES,0.0);
    tabEdge=vector<int>(MAXNODES);
    L = numsites;
    sameLTasPrev=vector<bool>(L,true);
    edge=vector<RecEdge*>(0);
    if(addedges) addEdgesFromFile(infile,0);
}

RecTree::RecTree(RecTree *intree,string newick,bool isFilename, bool forceages):Tree(newick,isFilename,forceages)
{
    tabNode=intree->tabNode;
    tabRec=intree->tabRec;
    tabSons=intree->tabSons;
    tabSonsDist=intree->tabSonsDist;
    tabFather=intree->tabFather;
    tabEdge=intree->tabEdge;
    age=intree->age;
    L = intree->L;
    edge=vector<RecEdge*>(intree->numRecEdge());
    sameLTasPrev=intree->sameLTasPrev;
    for(int i=0;i<intree->numRecEdge();i++) {
    	edge[i]= SlotAllocator<RecEdge>::GetSlotAllocator().Allocate();	// gets storage for a new RecEdge from a pool of storage
    	new (edge[i]) RecEdge(*(intree->getEdge(i)));	// uses in-place new to call copy constructor
    }
}

RecTree::RecTree(const RecTree& rt) :
	Tree(rt),
	tabSons(rt.tabSons),
	tabSonsDist(rt.tabSonsDist),
	tabNode(rt.tabNode),
	tabRec(rt.tabRec),
	tabFather(rt.tabFather),
	age(rt.age),
	tabEdge(rt.tabEdge),
	sameLTasPrev(rt.sameLTasPrev),
	L(rt.L)
{
    edge=vector<RecEdge*>(rt.numRecEdge());
    for(int i=0;i<rt.numRecEdge();i++) {
    	edge[i]= SlotAllocator<RecEdge>::GetSlotAllocator().Allocate();	// gets storage for a new RecEdge from a pool of storage
    	new (edge[i]) RecEdge(*(rt.getEdge(i)));	// uses in-place new to call copy constructor
    }
}

RecTree::RecTree(Data * data,double rho,double delta,vector<int> blocks):Tree(data)
{
    tabNode=vector<int>(n*2-1);
    tabRec=vector<int>(MAXNODES);
    tabSons=vector<vector<int> >(MAXNODES,vector<int>(2,0));
    tabSonsDist=vector<vector<double> >(MAXNODES,vector<double>(2,0.0));
    tabFather=vector<int>(MAXNODES,-1);
    age=vector<double>(MAXNODES,0.0);
    tabEdge=vector<int>(MAXNODES);
    L=data->getBlocks()->back();
    sameLTasPrev=vector<bool>(L,true);
    edge=vector<RecEdge*>(0);
    if(rho==0.0) return;else dropEdges(rho,delta,blocks);
}

RecTree::RecTree(int n,double rho,double delta,vector<int> blocks):Tree(n)
{
    tabNode=vector<int>(n*2-1);
    tabRec=vector<int>(MAXNODES);
    tabSons=vector<vector<int> >(MAXNODES,vector<int>(2,0));
    tabSonsDist=vector<vector<double> >(MAXNODES,vector<double>(2,0.0));
    tabFather=vector<int>(MAXNODES,-1);
    age=vector<double>(MAXNODES,0.0);
    tabEdge=vector<int>(MAXNODES);
    L=blocks.back();
    sameLTasPrev=vector<bool>(L,true);
    edge=vector<RecEdge*>(0);
    if(rho==0.0) return;else dropEdges(rho,delta,blocks);
}

void RecTree::assign(const RecTree& rt)
{
	Tree::assign(rt);
	//ensure everything is sized appropriately
	if(tabNode.size() != rt.tabNode.size())	tabNode.resize(rt.tabNode.size());
	if(tabRec.size() != rt.tabRec.size())	tabRec.resize(rt.tabRec.size());
	if(tabSons.size() != rt.tabSons.size())	tabSons.resize(rt.tabSons.size());
	if(tabSonsDist.size() != rt.tabSonsDist.size())	tabSonsDist.resize(rt.tabSonsDist.size());
	if(tabFather.size() != rt.tabFather.size())	tabFather.resize(rt.tabFather.size());
	if(age.size() != rt.age.size())	age.resize(rt.age.size());
	if(sameLTasPrev.size() != rt.sameLTasPrev.size())	sameLTasPrev.resize(rt.sameLTasPrev.size());
	unsigned int oldedgesize = edge.size();
	if(edge.size() != rt.edge.size())	edge.resize(rt.edge.size());
	//start copying data
	tabNode.assign(rt.tabNode.begin(), rt.tabNode.end());
	tabRec.assign(rt.tabRec.begin(), rt.tabRec.end());
	tabSons.assign(rt.tabSons.begin(), rt.tabSons.end());
	tabSonsDist.assign(rt.tabSonsDist.begin(), rt.tabSonsDist.end());
	tabFather.assign(rt.tabFather.begin(), rt.tabFather.end());
	age.assign(rt.age.begin(), rt.age.end());
	sameLTasPrev.assign(rt.sameLTasPrev.begin(), rt.sameLTasPrev.end());
	L=rt.L;
	unsigned int i=0;
	oldedgesize = oldedgesize < edge.size() ? oldedgesize : edge.size();
    for(;i<oldedgesize;i++) {
    	new (edge[i]) RecEdge(*(rt.getEdge(i)));	// uses in-place new to call copy constructor
    }
    for(;i<edge.size();i++) {
    	edge[i]= SlotAllocator<RecEdge>::GetSlotAllocator().Allocate();	// gets storage for a new RecEdge from a pool of storage
    	new (edge[i]) RecEdge(*(rt.getEdge(i)));	// uses in-place new to call copy constructor
    }
}


void RecTree::addEdgesFromFile(WargXml *infile,int siteoffset)
{
	long numedges=0;
	streampos cp=infile->tellg();
	while(1) {
		string res=infile->getLine();
		if(res.find("<recedge>")!=string::npos){
			try{
			int pos=addRecEdge(res,siteoffset);
			if(pos>=0) numedges++;
			else if(pos==-1)cerr<<"WARNING: failed to add an observed edge!"<<endl;
			}catch(char * x){cerr<<"Error in reading XML file:"<<x<<endl;exit(0);}
		}
		if (infile->eof() || res.find("</Iteration>")!=string::npos) break;
	}
	infile->clear();
	infile->seekg(cp);
	
	dlog(1)<<"Added "<<numedges<<" Recombination edge(s) from output file"<<endl;

}


void RecTree::dropEdges(double rho,double delta,vector<int> blocks)
{
    double tfrom,tto;
    unsigned int gstart,gend,edgefrom,edgeto;
    vector<int> nl;
    for(int i=0;i<n;i++)
        nl.push_back(nodes[i]->getFather()->getId());// don't actually use fatherId at the mo
    for(int i=n;i<n+n-2;i++)
        nl.push_back(-1);// -1 means not alive now
    int k=n;
    double time=0.0;
    while (k>1)
    {
#if defined DEBUG
        for(int i=0;i<n+n-2;i++)
            cout<<"," <<nl[i];
        cout << ": k="<<k<<" time="<<time<<endl;
#endif

        time-=2.0/(rho * k)*log(gsl_rng_uniform(rng));
        if(time >  nodes[n+n-k]->getAge())
        {// no recedge occurred before next coalescence
            time=nodes[n+n-k]->getAge();
#if defined DEBUG

            cout << "COALESCENCE: time="<<time<<endl;
#endif

            if(k>2)
            {
                nl[nodes[n+n-k]->getLeft()->getId()]=-1;
                nl[nodes[n+n-k]->getRight()->getId()]=-1;
                nl[n+n-k]=nodes[n+n-k]->getFather()->getId();
            }
            k--;
        }
        else
        {// recombination event
            tto=time;
            tfrom=time;//note:updated in getEdgeCoal
            edgeto= randomActiveNode(nl);
            edgefrom=getEdgeCoal(nl,k,&tfrom);
            setBlock(&gstart,&gend,delta,&blocks);
#if defined DEBUG

            cout <<"tfrom="<<tfrom<<" tto=" <<tto <<" gstart="<<gstart<<" gend="<< gend <<" edgefrom="<<edgefrom<<" edgeto="<<edgeto<<endl;
#endif

            if (tfrom<tto)
            {
                cerr<<"Error in simulation"<<endl;
                throw "Edge arrives before departs";
            }
            if(addRecEdge(tfrom-nodes[edgefrom]->getAge(),tto-nodes[edgeto]->getAge(),gstart,gend,edgefrom,edgeto)<0) throw("RecTree: Can't create tree");
        }
    }
#if defined DEBUG
    testTree();
#endif
}

void RecTree::setBlock(unsigned int* gstart,unsigned int* gend,double delta,vector<int>* blocks)
{
    unsigned int b=blocks->size()-1, L=blocks->back(), blockin=0;
    double tmp=gsl_rng_uniform(rng);
    double pofblock=0;
    // prob of block \propto (delta-1+length of block)
    while(blockin < b-1 && tmp> (pofblock+= double(delta-1+(blocks->at(blockin+1)-blocks->at(blockin)))/(b*(delta-1)+L)))
    {
        blockin++;
    }// determine which block is affected
    // prob starting within block|block i = (length of block-1)/(delta+length of block-1)
    if(gsl_rng_uniform(rng)<(blocks->at(blockin+1)-blocks->at(blockin)-1)/(delta+(blocks->at(blockin+1)-blocks->at(blockin)-1)))
    {// starts within block
        *gstart=blocks->at(blockin)+1+(int)floor(gsl_rng_uniform(rng)*(blocks->at(blockin+1)-blocks->at(blockin)-1));
    }
    else
    {
        *gstart=blocks->at(blockin);
    }
    if((*gend=*gstart+gsl_ran_geometric(rng,1.0/delta))>= (unsigned int)blocks->at(blockin+1))
        *gend=blocks->at(blockin+1);// sample from the geometric distribution with minimum size 1 (so minimum delta g=0)
}

int RecTree::randomActiveNode(vector<int> nodelist)
{
    vector<int> livenodes;
    for(unsigned int i=0;i<nodelist.size();i++)
    {
        if(nodelist[i]>=0)
            livenodes.push_back(i);
    }
    return livenodes[gsl_rng_uniform_int(rng,livenodes.size())];
}

int RecTree::getEdgeCoal(vector<int> nl,int k,double* time)
{
    for(;k>1;k--)
    {
        *time-=log(gsl_rng_uniform(rng))/k; // any given edge coaelesces at rate 1 WITH ALL OTHER EXTANT EDGES
        if(*time<nodes[n+n-k]->getAge())
        {// coalescence of this edge occurs
            return randomActiveNode(nl);
        }
        else
        {// coalescence within the tree occurs
            *time=nodes[n+n-k]->getAge();
            nl[nodes[n+n-k]->getLeft()->getId()]=-1;
            nl[nodes[n+n-k]->getRight()->getId()]=-1;
            if(k>2)
                nl[n+n-k]=nodes[n+n-k]->getFather()->getId();
        }
    }
    *time-=log(gsl_rng_uniform(rng));
    return n+n-k-1; // when only one type, recedge will coalesce with it
}

int RecTree::getEdgeCoal(double* time)
{
    vector<int> nl(n+n-2);
    for(int i=0;i<n;i++)
        nl[i]=nodes[i]->getFather()->getId();
    for(int i=n;i<n+n-2;i++)
        nl[i]=-1;
    int k=n;
    while (k>1)
    {
        if (*time >  nodes[n+n-k]->getAge())
        {
            nl[nodes[n+n-k]->getLeft()->getId()]=-1;
            nl[nodes[n+n-k]->getRight()->getId()]=-1;
            nl[n+n-k]=nodes[n+n-k]->getFather()->getId();
        }
        else
            break;
        k--;
    }
    return getEdgeCoal(nl,k,time);
}

RecTree::~RecTree()
{
    for (unsigned int i=0;i<edge.size();i++)
        if(edge[i]!=NULL){
        	SlotAllocator<RecEdge>::GetSlotAllocator().Free(edge[i]);
        }
}

int RecTree::addRecEdge(double tfrom,double tto,unsigned int gstart,unsigned int gend,int edgefrom,int edgeto)
{
    // Need to add some extra error handling routines here
///*** NOTE: There has been an awful hack implemented in here, ignoring edges that would otherwise be illegal.  IF YOU GET WARNINGS OR ERRORS, AND YOU ARE NOT READING IN A TRUNCATED RECEDGES FROM A FILE, YOU SHOULD BE VERY WORRIED!  If you are reading in a recedge with low precision, then it can be ignored, which is bad but usually not terrible.
    int N=getN();
try{
    if(edgeto>2 * N - 3)
    {
        cerr << "Edge seems to arrive above the root!" <<endl;
        throw "Edge root arrival violation";
    }
    if(edgefrom>2 * N - 2)
    {
        cerr << "Edge seems to depart outside the tree!" <<endl;
        throw "Edge tree depart violation";
    }
    if(edgeto  <2*N-2&&tto  +nodes[edgeto  ]->getAge()>nodes[edgeto  ]->getFather()->getAge())
    {
        cerr <<std::setprecision(64)<< "Edge to age:"<<tto+nodes[edgeto  ]->getAge()<<" appears to be greater than that nodes parent age:"<<nodes[edgeto]->getFather()->getAge()<<". Edge details: tfrom="<<tfrom<<" tto="<<tto<<" efrom="<<edgefrom<<" eto="<<edgeto<<" start="<<gstart<<" end="<<gend<<endl;
	if(tto  +nodes[edgeto  ]->getAge()-nodes[edgeto  ]->getFather()->getAge()>0.0001)        throw "Edge age violation";// throw error if difference is not small
	cerr<<"Error in addRecEdge: Not adding edge!"<<endl<<endl;
	return(-1);// or just ignore edge otherwise;
    }
    if(edgefrom<2*N-2&&tfrom+nodes[edgefrom]->getAge()>nodes[edgefrom]->getFather()->getAge())
    {
	// allow for some loss of precision
	if((tfrom+nodes[edgefrom]->getAge())-nodes[edgefrom]->getFather()->getAge()>0.0001)
	{
	        cerr <<std::setprecision(64)<< "Edge from node "<<edgefrom<<" of age "<<nodes[edgefrom]->getAge()<< " has age:"<<tfrom+nodes[edgefrom]->getAge()<<" which is greater than that nodes parent age:"<<nodes[edgefrom]->getFather()->getAge()<<"! (from "<<edgefrom<<")"<<endl;
	        throw "Edge age violation";
	}else if((tfrom+nodes[edgefrom]->getAge())-nodes[edgefrom]->getFather()->getAge() > 0){
		//tfrom=(nodes[edgefrom]->getFather()->getAge()-nodes[edgefrom]->getAge())*(1.0-0.001 * gsl_rng_uniform(rng));
		cerr<<"Error in addRecEdge: Not adding edge!"<<endl<<endl;
		return(-1);// or just ignore edge otherwise;
	}
        if(gstart==gend)
        {
            cerr<<"gstart==gend"<<endl;
            throw "Empty recombination error";
        }
    }
}catch(char * x){
	cerr<<"Error in addRecEdge: Not adding edge!"<<endl<<x<<endl;
	return(-1);
}
    unsigned int i=0;
    while (i<edge.size() && edge[i]->getTimeFrom()+nodes[edge[i]->getEdgeFrom()]->getAge()<tfrom+nodes[edgefrom]->getAge())
        i++;

    RecEdge* re = SlotAllocator<RecEdge>::GetSlotAllocator().Allocate();	// gets storage for a new RecEdge from a pool of storage
	new (re) RecEdge(tfrom,tto,gstart,gend,edgefrom,edgeto);	// uses in-place new to call constructor

    edge.insert(edge.begin()+i,re);
    sameLTasPrev[gstart]=false;
    if (gend<L)
        sameLTasPrev[gend]=false;
    return i;
}

int RecTree::addRecEdge(std::string res,int sitesoffset)
{
    double tfrom,tto;
    int gstart,gend;
    int edgefrom,edgeto;
    size_t f1,f2;
    f1=res.find("<start>")+7;
    f2=res.find("</start>",f1);
    gstart=atoi(res.substr(f1,f2-f1).c_str())+sitesoffset;
    f1=res.find("<end>",f2)+5;
    f2=res.find("</end>",f1);
    gend=atoi(res.substr(f1,f2-f1).c_str())+sitesoffset;
    f1=res.find("<efrom>",f2)+7;
    f2=res.find("</efrom>",f1);
    edgefrom=atoi(res.substr(f1,f2-f1).c_str());
    f1=res.find("<eto>",f2)+5;
    f2=res.find("</eto>",f1);
    edgeto=atoi(res.substr(f1,f2-f1).c_str());
    f1=res.find("<afrom>",f2)+7;
    f2=res.find("</afrom>",f1);
    tfrom=atof(res.substr(f1,f2-f1).c_str());
    f1=res.find("<ato>",f2)+5;
    f2=res.find("</ato>",f1);
    tto=atof(res.substr(f1,f2-f1).c_str());
    if(gend>=(int)L || gstart<0) return(-2);
    return addRecEdge(tfrom,tto,gstart,gend,edgefrom,edgeto);
}


int RecTree::addRecEdge(unsigned int gstart,unsigned int gend,int edgefrom,int edgeto)
{
    double tfrom=gsl_rng_uniform(rng) * getNode(edgefrom)->getDist();
    double tto=gsl_rng_uniform(rng) * getNode(edgeto)->getDist();
    return addRecEdge(tfrom,tto,gstart,gend,edgefrom,edgeto);
}

void RecTree::remRecEdge(int which)
{
    unsigned int start=edge[which]->getStart();
    unsigned int end=edge[which]->getEnd();
    if(edge[which]!=NULL) {
    	SlotAllocator<RecEdge>::GetSlotAllocator().Free(edge[which]);
    }
    edge.erase(edge.begin()+which);
    if (start>0)
        calcSameLTasPrev(start);
    if (end<L)
        calcSameLTasPrev(end);
}


/** Fast exponential taken from Cawley 2000, Neural Computation */
#define EXPA (1048576/M_LN2)
#define EXPC 60801
inline double exponential(double y)
{
   union
   {
       double d;
#ifdef LITTLE_ENDIAN
       struct { int j, i; } n;
#else
       struct { int i, j; } n;
#endif
   }
   eco;
   eco.n.i = (int)(EXPA*(y)) + (1072693248 - EXPC);
   eco.n.j = 0;
   return eco.d;
}
/** end taken from Cawley 2000 */


void RecTree::makeLocalTree(unsigned int site,double thetaPerSite)
{
    //Sort nodes and departure points in increasing order of age
    int indTabEdge=0;
    int aff=affecting(site);
    for (int i=0;i<n;i++)
        tabNode[i]=i;
    unsigned int indNode=n;
    unsigned int indEdge=0;
    bool isNode=true;
    for (int j=0;j<n-1+aff;j++)
    {
        while (indEdge<edge.size() && edge[indEdge]->affectsSite(site)==false)
            indEdge++;
        if (indNode>=nodes.size())
        {
            isNode=false;
            indEdge++;
            goto next;
        };
        if (indEdge>=edge.size())
        {
            isNode=true;
            indNode++;
            goto next;
        }
        if (nodes[indNode]->getAge()<getEdgeTimeAbsFrom(indEdge))
        {
            isNode=true;
            indNode++;
        }
        else
        {
            isNode=false;
            indEdge++;
        }
next:
        if (isNode)
        {
            tabNode[indNode-1]=j+n;
            age[j+n]=nodes[indNode-1]->getAge();
        }
        else
        {
            tabRec [indEdge-1]=j+n;
            age[j+n]=getEdgeTimeAbsFrom(indEdge-1);
            tabEdge[indTabEdge++]=indEdge-1;
        }
    }
    //Find children of each node
    for (int i=n;i<2*n-1+aff;i++)
        for (int isLeft=0;isLeft<=1;isLeft++)
        {
            //Find out which edge "ed" to go down from and starting from age "ageCurrent"
            int ind,ii;
            double ageCurrent;
            int ed;
            if (i<2*n-1)
            {
                isNode=true;
                ind=i;
                ii=tabNode[ind];
                ageCurrent=nodes[ind]->getAge();
                if (isLeft)
                {
                    ed=nodes[ind]->getLeft ()->getId();
                }
                else
                {
                    ed=nodes[ind]->getRight()->getId();
                }
            }
            else
            {
                isNode=false;
                ind=tabEdge[i-2*n+1];
                ii=tabRec[ind];
                if (isLeft)
                {
                    ed=edge[ind]->getEdgeFrom();
                    ageCurrent=getEdgeTimeAbsFrom(ind);
                }
                else
                {
                    ed=edge[ind]->getEdgeTo();
                    ageCurrent=getEdgeTimeAbsTo(ind);
                }
            }

            //Go down the edge to find the son (if any)
            int cur=-1;
            double curAge=0;
            bool from=false;
            for (int k=0;k<aff;k++)
            {
                if (edge[tabEdge[k]]->getEdgeFrom()==ed && isless(getEdgeTimeAbsFrom(tabEdge[k]),ageCurrent) && (cur==-1||isless(curAge,getEdgeTimeAbsFrom(tabEdge[k]))))
                {
                    cur=tabEdge[k];
                    curAge=getEdgeTimeAbsFrom(tabEdge[k]);
                    from=true;
                };
                if (edge[tabEdge[k]]->getEdgeTo  ()==ed && isless(getEdgeTimeAbsTo  (tabEdge[k]),ageCurrent) && (cur==-1||isless(curAge,getEdgeTimeAbsTo  (tabEdge[k]))))
                {
                    cur=tabEdge[k];
                    curAge=getEdgeTimeAbsTo  (tabEdge[k]);
                    from=false;
                };
            }
            if (cur==-1)
            {
                isNode=true;
                ind=ed;
            }
            else if (from)
            {
                isNode=false;
                ind=cur;
            }
            else
            {
                tabSons[ii][isLeft]=-1;
                continue;
            }
            if (isNode)
                tabSons[ii][isLeft]=tabNode[ind];
            else
                tabSons[ii][isLeft]=tabRec[ind];
            tabFather[tabSons[ii][isLeft]]=ii;
        }
    //Remove unnecessary nodes
    int cur=n;
    for (int i=n;i<2*n-1+aff;i++)
    {
        int f=tabFather[i];
        if (tabSons[i][0]==-1 && tabSons[i][1]==-1)
        {
            if (tabSons[f][0]==i)
                tabSons[f][0]=-1;
            else
                tabSons[f][1]=-1;
            continue;
        }
        if (tabSons[i][0]==-1 && f>=0)
        {
            if (tabSons[f][0]==i)
                tabSons[f][0]=tabSons[i][1];
            else
                tabSons[f][1]=tabSons[i][1];
            tabSons[i][1]=-1;
            continue;
        };
        if (tabSons[i][1]==-1 && f>=0)
        {
            if (tabSons[f][0]==i)
                tabSons[f][0]=tabSons[i][0];
            else
                tabSons[f][1]=tabSons[i][0];
            tabSons[i][0]=-1;
            continue;
        };
        if (tabSons[i][0]>=0 && tabSons[i][1]>=0)
        {
            tabSons[cur][0]=tabSons[i][0];
            tabSons[cur][1]=tabSons[i][1];
            age[cur]=age[i];
            if (f>=0)
            {
                if (tabSons[f][0]==i)
                    tabSons[f][0]=cur;
                else
                    tabSons[f][1]=cur;
            }

//            tabSonsDist[cur][0]=exponential(-(2.0/3.0)*(age[cur]-age[tabSons[cur][0]])*thetaPerSite);
//            tabSonsDist[cur][1]=exponential(-(2.0/3.0)*(age[cur]-age[tabSons[cur][1]])*thetaPerSite);
            tabSonsDist[cur][0]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][0]])*thetaPerSite);
            tabSonsDist[cur][1]=exp(-(2.0/3.0)*(age[cur]-age[tabSons[cur][1]])*thetaPerSite);
            if(tabSonsDist[cur][0]>1||tabSonsDist[cur][1]>1)
            {
                cerr<<std::setprecision(64)<<"makeLocalTree:negative sonsdist! "<<age[cur]<<"-"<<age[tabSons[cur][0]]<<" and "<<age[cur]<<"-"<<age[tabSons[cur][1]]<<endl;
                throw("Negative sons dist");
            }
            cur++;
        }
    }
}


void RecTree::makeLocalTreeKeepingFathers(unsigned int site)
{
    //Sort nodes and departure points in increasing order of age
    int indTabEdge=0;
    int aff=affecting(site);
    for (int i=0;i<n;i++)
        tabNode[i]=i;
    unsigned int indNode=n;
    unsigned int indEdge=0;
    bool isNode=true;
    for (int j=0;j<n-1+aff;j++)
    {
        while (indEdge<edge.size() && edge[indEdge]->affectsSite(site)==false)
            indEdge++;
        if (indNode>=nodes.size())
        {
            isNode=false;
            indEdge++;
            goto next;
        };
        if (indEdge>=edge.size())
        {
            isNode=true;
            indNode++;
            goto next;
        }
        if (nodes[indNode]->getAge()<getEdgeTimeAbsFrom(indEdge))
        {
            isNode=true;
            indNode++;
        }
        else
        {
            isNode=false;
            indEdge++;
        }
next:
        if (isNode)
        {
            tabNode[indNode-1]=j+n;
            age[j+n]=nodes[indNode-1]->getAge();
        }
        else
        {
            tabRec [indEdge-1]=j+n;
            age[j+n]=getEdgeTimeAbsFrom(indEdge-1);
            tabEdge[indTabEdge++]=indEdge-1;
        }
    }
    //Find children of each node
    for (int i=n;i<2*n-1+aff;i++){
        for (int isLeft=0;isLeft<=1;isLeft++)
        {
            //Find out which edge "ed" to go down from and starting from age "ageCurrent"
            int ind,ii;
            double ageCurrent;
            int ed;
            if (i<2*n-1)
            {
                isNode=true;
                ind=i;
                ii=tabNode[ind];
                ageCurrent=nodes[ind]->getAge();
                if (isLeft)
                {
                    ed=nodes[ind]->getLeft ()->getId();
                }
                else
                {
                    ed=nodes[ind]->getRight()->getId();
                }
            }
            else
            {
                isNode=false;
                ind=tabEdge[i-2*n+1];
                ii=tabRec[ind];
                if (isLeft)
                {
                    ed=edge[ind]->getEdgeFrom();
                    ageCurrent=getEdgeTimeAbsFrom(ind);
                }
                else
                {
                    ed=edge[ind]->getEdgeTo();
                    ageCurrent=getEdgeTimeAbsTo(ind);
                }
            }

            //Go down the edge to find the son (if any)
            int cur=-1;
            double curAge=0;
            bool from=false;
            for (int k=0;k<aff;k++)
            {
                if (edge[tabEdge[k]]->getEdgeFrom()==ed && isless(getEdgeTimeAbsFrom(tabEdge[k]),ageCurrent) && (cur==-1||isless(curAge,getEdgeTimeAbsFrom(tabEdge[k]))))
                {
                    cur=tabEdge[k];
                    curAge=getEdgeTimeAbsFrom(tabEdge[k]);
                    from=true;
                };
                if (edge[tabEdge[k]]->getEdgeTo  ()==ed && isless(getEdgeTimeAbsTo  (tabEdge[k]),ageCurrent) && (cur==-1||isless(curAge,getEdgeTimeAbsTo  (tabEdge[k]))))
                {
                    cur=tabEdge[k];
                    curAge=getEdgeTimeAbsTo  (tabEdge[k]);
                    from=false;
                };
            }
            if (cur==-1)
            {
                isNode=true;
                ind=ed;
            }
            else if (from)
            {
                isNode=false;
                ind=cur;
            }
            else
            {
                tabSons[ii][isLeft]=-1;
                continue;
            }
            if (isNode)
                tabSons[ii][isLeft]=tabNode[ind];
            else
                tabSons[ii][isLeft]=tabRec[ind];
            tabFather[tabSons[ii][isLeft]]=ii;
        }
    }
}


void RecTree::testEdges() const
{
    try
    {
        for(int i=0;i<numRecEdge();i++)
        {
            if(i>0)
            {
                if(getEdgeTimeAbsFrom(i-1)>getEdgeTimeAbsFrom(i))
                {
                    cerr<<std::setprecision(64)<<"Edge "<<i-1<<" is from time "<<getEdgeTimeAbsFrom(i-1)<<" which is greater than the next edge from time "<<getEdgeTimeAbsFrom(i)<< endl;
                    throw "Edge ordering incorrect";
                }
            }
            if(getEdgeTimeAbsFrom(i)<getEdgeTimeAbsTo(i))
            {
                cerr<<std::setprecision(64)<<"Edge "<<i<<" is from time "<<getEdgeTimeAbsFrom(i)<<" which is less than the time it goes to at "<<getEdgeTimeAbsTo(i)<< endl;
                throw "Edge direction violated";
            }
            if(getNode(getRecEdge(i)->getEdgeFrom())->getFather()!=NULL)
            {
                if(getNode(getRecEdge(i)->getEdgeFrom())->getFather()->getAge()< getEdgeTimeAbsFrom(i))
                {
                    cerr<<std::setprecision(64)<<"Edge "<<i<<" is from "<<getEdgeTimeAbsFrom(i)<<" which is greater than nodes fathers age "<< getNode(getRecEdge(i)->getEdgeFrom())->getFather()->getAge()<< endl;
                    throw "Illegal edge";
                }
            }
            if(getNode(getRecEdge(i)->getEdgeFrom())->getAge()>getEdgeTimeAbsFrom(i))
            {
                cerr<<std::setprecision(64)<<"Edge "<<i<<" is from "<<getEdgeTimeAbsFrom(i)<<" which is less than nodes age "<< getNode(getRecEdge(i)->getEdgeFrom())->getAge()<< endl;
                throw "Illegal edge";
            }
            if(getNode(getRecEdge(i)->getEdgeTo())->getFather()!=NULL)
            {
                if(getNode(getRecEdge(i)->getEdgeTo())->getFather()->getAge()< getEdgeTimeAbsTo(i))
                {
                    cerr<<std::setprecision(64)<<"Edge "<<i<<" is to "<<getEdgeTimeAbsTo(i)<<" which is greater than nodes fathers age "<< getNode(getRecEdge(i)->getEdgeTo())->getFather()->getAge()<< endl;
                    throw "Illegal edge";
                }
            }
            if(getNode(getRecEdge(i)->getEdgeTo())->getAge()>getEdgeTimeAbsTo(i))
            {
                cerr<<std::setprecision(64)<<"Edge "<<i<<" is to "<<getEdgeTimeAbsTo(i)<<" which is less than nodes age "<< getNode(getRecEdge(i)->getEdgeTo())->getAge()<< endl;
                throw "Illegal edge";
            }
        }
        dlog(1)<<"Testing recedges complete"<<endl;
        //throw "OK";// this is to always get the debug info
    }
    catch (char * x)
    {
        // Debugging info
  /*         for(int i=0;i<numRecEdge();i++)
        {
         cout<<" Edge "<<i<<" from Node "<< getRecEdge(i)->getEdgeFrom()<< "("<< getNode(getRecEdge(i)->getEdgeFrom())->getAge() <<":";
            if(getNode(getRecEdge(i)->getEdgeFrom())->getFather()!=NULL )
                cout << getNode(getRecEdge(i)->getEdgeFrom())->getFather()->getAge();
            cout<<")"<<" at time " << getEdgeTimeAbsFrom(i)<<" to "<< getRecEdge(i)->getEdgeTo()<< "("<< getNode(getRecEdge(i)->getEdgeTo())->getAge() <<":";
            if(getNode(getRecEdge(i)->getEdgeTo())->getFather()!=NULL)
                cout<<getNode(getRecEdge(i)->getEdgeTo())->getFather()->getAge();
            cout<<") at time "<< getEdgeTimeAbsTo(i)<<endl;
        }*/
        if(strcmp(x,"OK"))
            exit(1);
    }
}

void RecTree::testTree()
{
	dlog(1)<<"Testing tree"<<endl;
    testEdges();
    dlog(1)<<"Testing local trees"<<endl;
    for(unsigned int site =0;site<L;site++)
    {
        makeLocalTree(site,1.0);
        for (int i=n;i<2*n-1;i++)
            for (int side=0;side<2;side++)
            {
                double dist;
                int j=getSon(i,site,side,&dist);
                if (j>=0&&j>=i)
                {
                    cerr<<std::setprecision(64)<<"Error: "<<i<<" (age "<<nodes[i]->getAge()<<") has son "<<j<<" (age "<<nodes[j]->getAge()<<")"<<endl;
                    for (unsigned int i=0;i<edge.size();i++)
                        if (edge[i]->affectsSite(site))
                            cout<<std::setprecision(64) <<"Edge "<<i<<" goes from "<<edge[i]->getEdgeFrom()<<":"<<getEdgeTimeAbsFrom(i)<<" to "<<edge[i]->getEdgeTo()<<":"<<getEdgeTimeAbsTo(i)<<endl;
                    throw "testtree: Impossible son";
                };
                if (j>=0&&dist<=0.0)
                {
                    cerr<<std::setprecision(64)<<"Error with dist="<<dist<<endl;
                    throw "testtree: negative distance";
                }
            }
    }
    dlog(1)<<"Test complete."<<endl;
}

double RecTree::priorEdge(int e,Param * param) const
{
    double timeFrom=getEdgeTimeAbsFrom(e);
    double timeTo  =getEdgeTimeAbsTo  (e);
    double L=0.0;
    for (unsigned int i=n;i<nodes.size();i++)
    {
        if (timeFrom<nodes[i-1]->getAge() || timeTo>nodes[i]->getAge())
            continue;
        L+=(2*n-i)*(min(timeFrom,nodes[i]->getAge())-max(timeTo,nodes[i-1]->getAge()));
        if (timeFrom<nodes[i]->getAge())
            break;
    }
    if (timeFrom>nodes.back()->getAge())
        L+=timeFrom-nodes.back()->getAge();
    double ret=-L;
    int x=edge[e]->getStart();
    int y=edge[e]->getEnd();
    bool xatborder=false,yatborder=false;
    vector<int>*blocks=param->getData()->getBlocks();
    for (unsigned int i=0;i<blocks->size();i++) {
        if (blocks->at(i)==x) xatborder=true;
        if (blocks->at(i)==y) yatborder=true;
    }
    int b=blocks->size()-1;
    int length=blocks->back();
    double delta=param->getDelta();
    if (xatborder) ret+=log(delta/(b*delta+length-b)); else ret+=log(1.0/(b*delta+length-b));
    if (yatborder) ret+=(y-x+1)*log(1.0-1.0/delta); else ret+=(y-x)*log(1.0-1.0/delta)-log(delta);
    return(ret);
}

double RecTree::priorEdge(double tFrom,double tTo) const
{
    double L=0.0;
    if (tFrom<=tTo) return 0.0;
    for (unsigned int i=n;i<nodes.size();i++)
    {
        if (tFrom<nodes[i-1]->getAge() || tTo>nodes[i]->getAge())
            continue;

L+=(2*n-i)*(min(tFrom,nodes[i]->getAge())-max(tTo,nodes[i-1]->getAge()));
        if (tFrom<nodes[i]->getAge())
            break;
    }
    if (tFrom>nodes.back()->getAge())
        L+=tFrom-nodes.back()->getAge();
    return exp(-L)/ttotal;
}

double RecTree::prior(Param * param) const
{
    //Prior probability of clonal genealogy
    double ret=((Tree*)this)->prior();
    //Number of edges
    if (param->getRho()>0.0) {
    double T=getTTotal();
    ret+=-param->getRho()*T*0.5+edge.size()*log(param->getRho()*0.5);} else if (edge.size()>0) ret=-INFINITY;
    //Location of the edges
    for (unsigned int e=0;e<edge.size();e++)
        ret+=priorEdge(e,param);
    return ret;
}

void RecTree::calcSameLTasPrev(unsigned int site)
{
    for (unsigned int i=0;i<edge.size();i++)
        if (edge[i]->getStart()==site || edge[i]->getEnd()==site)
        {
            sameLTasPrev[site]=false;
            return;
        };
    sameLTasPrev[site]=true;
}

vector<int> RecTree::alive(double t) const
{
    vector<int> v;
    for (unsigned int i=0;i<nodes.size()-1;i++)
        if (nodes[i]->getAge()<=t && nodes[i]->getFather()->getAge()>=t)
            v.push_back(i);
    if (nodes.back()->getAge()<=t)
        v.push_back(nodes.size()-1);
    return v;
}

bool RecTree::isThere(int start,int end,int e,double time1,double time2) const
{
    if (time1>time2)
        swap(time1,time2);
    for (unsigned int i=0;i<edge.size();i++)
    {
        if (start>=(int)edge[i]->getEnd())
            continue;
        if ((int)edge[i]->getStart()>=end)
            continue;
        if (edge[i]->getEdgeFrom()==e && edge[i]->getTimeFrom()>=time1 && edge[i]->getTimeFrom()<=time2)
            return true;
        if (edge[i]->getEdgeTo  ()==e && edge[i]->getTimeTo  ()>=time1 && edge[i]->getTimeTo  ()<=time2)
            return true;
    }
    return false;
}

void RecTree::setStart(int e,int start)
{
    if(e>=(int)edge.size()||e<0) return;
    if(start<0 || start>(int)sameLTasPrev.size()) return;
    int old=edge[e]->gstart;
    edge[e]->gstart=start;
    if (old>0)
        calcSameLTasPrev(old);
    if (start>0)
        sameLTasPrev[start]=0;
}

void RecTree::setEnd(int e,int end)
{
    if(e>=(int)edge.size()||e<0) return;
    if(end<0 || end>(int)sameLTasPrev.size()) return;

    int old=edge[e]->gend;
    edge[e]->gend=end;
    if (old<(int)L)
        calcSameLTasPrev(old);
    if (end<(int)L)
        sameLTasPrev[end]=0;
}

void RecTree::scaleEdges(int which, double dist){
	// for every recedge going to which, move it proportional to dist
	double reldist;
	if(getNode(which)->getDist()>0)reldist=(getNode(which)->getDist() + dist)/getNode(which)->getDist();
	else reldist=0;
	for(int e=0;e< numRecEdge();e++) {
		if(edge[e]->getEdgeTo()==which)  updateEdgeTimes(e,reldist);
	}
}

void RecTree::fixToTimes(int which, double dist){
	// for every recedge going to which, move it proportional to dist
// NOTE: which *can* be the root, because the tree may not be correctly ordered at this point.
	for(int e=0;e< numRecEdge();e++) {
		if(edge[e]->getEdgeTo()==which)  {
			edge[e]->setTimeTo(edge[e]->getTimeTo()+dist);
			edge[e]->setTimeFrom(edge[e]->getTimeFrom()+dist);
		}
	}
}

void RecTree::fixFromTimes(int affedge,double dist)
{
	for(unsigned int i=0;i<edge.size();i++) {
		if(edge[i]->getEdgeTo()==affedge){
			edge[i]->setTimeFrom(edge[i]->getTimeFrom()+dist);
		}
		if(edge[i]->getEdgeFrom()==affedge){
			edge[i]->setTimeFrom(edge[i]->getTimeFrom()-dist);
		}
	}
}

void RecTree::changeAge(int which,double dist)
{
	fixFromTimes(which,dist);
	getNode(which)->changeAge(dist);
}

void RecTree::swapEdgeTo(int a, int b)
{
	double dist;
	for(unsigned int i=0;i<edge.size();i++) {
		if(edge[i]->getEdgeTo()==a)
		{
			dist=getEdgeTimeAbsFrom(i)-getEdgeTimeAbsTo(i);
			edge[i]->setTimeTo(getEdgeTimeAbsTo(i)-getNode(b)->getAge(),b);
			edge[i]->setTimeFrom(getEdgeTimeAbsTo(i)+dist);
		}else if(edge[i]->getEdgeTo()==b) {
			dist=getEdgeTimeAbsFrom(i)-getEdgeTimeAbsTo(i);
			edge[i]->setTimeTo(getEdgeTimeAbsTo(i)-getNode(a)->getAge(),a);
			edge[i]->setTimeFrom(getEdgeTimeAbsTo(i)+dist);
		}
	}
}

void RecTree::chooseNodePath(int which, double dist, vector<int> *listLR, vector<int> *affedges)
{
	int  newfather,notnewfather,oldfather;
	vector<int> nnflist;
	vector<double> test;
	for(int i=0;i<(int)nodes.size();i++) test.push_back(nodes[i]->getAge());
	changeAge(which,dist);//change the age of a node (maintaining recedges)
	affedges->push_back(which);
	affedges->push_back(getNode(which)->getLeft()->getId());
	affedges->push_back(getNode(which)->getRight()->getId());
	int revnode=getOldestReversedNode();
	while(revnode>=0){// construct the list of affedges
		if(gsl_rng_uniform(rng)>0.5) {
			newfather=getNode(revnode)->getLeft()->getId();
			notnewfather=getNode(revnode)->getRight()->getId();
			listLR->push_back(1);
		}else{
			notnewfather=getNode(revnode)->getLeft()->getId();
			newfather=getNode(revnode)->getRight()->getId();
			listLR->push_back(0);
		}
		oldfather=getNode(revnode)->getFather()->getId();
		swapFather(notnewfather,revnode);
		swapFather(revnode,oldfather);
		affedges->push_back(oldfather);
		nnflist.push_back(notnewfather);
		revnode=getOldestReversedNode();
	}
	// put things back as they were
	for(int i=listLR->size()-1;i>=0;i--) {
	  revnode=getNode(affedges->at(3+i))->getId();
	  oldfather=getNode(revnode)->getFather()->getId();
	  notnewfather=nnflist[i];
	  if(getNode(revnode)->getRight()->getId()==notnewfather) newfather=getNode(revnode)->getLeft()->getId();
	  else if(getNode(revnode)->getLeft()->getId()==notnewfather) newfather=getNode(revnode)->getRight()->getId();
	  else throw("Incompatible changes!");
	  if(notnewfather==newfather)throw("Notnewfather == newfather!");
	  swapFather(notnewfather,revnode);
	  swapFather(revnode,oldfather);
	}
	changeAge(which,-dist);
	for(int i=0;i<(int)nodes.size();i++) if(fabs(test[i]-nodes[i]->getAge())>0.0001) throw("Caught a bug!");
}

int RecTree::moveNodeTime(int which, int *whichto,double dist,int oldwhich, vector<int> *listLR,  vector<int> *affedges)
{
	// scale the recedge events
	// scale (first, so that the branch lengths are simple to work with)
	scaleEdges(which,-dist);
	scaleEdges(getNode(which)->getLeft()->getId(),dist);
	scaleEdges(getNode(which)->getRight()->getId(),dist);
	changeAge(which,dist);//change the age of a node (maintaining recedges)

	int  newfather,notnewfather;
	int revnode=getOldestReversedNode();

	for(int i=0;i<(int)listLR->size();i++) {
		if(revnode<0) throw("RecTree: Incompatible list error");
		// choose which of the daughters gets promoted to the parent
		if(listLR->at(i)==1) {
			newfather=getNode(revnode)->getLeft()->getId();
			notnewfather=getNode(revnode)->getRight()->getId();
		}else{
			notnewfather=getNode(revnode)->getLeft()->getId();
			newfather=getNode(revnode)->getRight()->getId();
		}
		// Second round of scaling: revnodes father, and the edge that is not the new father
		double newdist=getNode(revnode)->getAge()-getNode(revnode)->getFather()->getAge();
		int oldfather=getNode(revnode)->getFather()->getId();
		scaleEdges(oldfather,-newdist);
		scaleEdges(notnewfather,-newdist);
		// we keep absolute time constant in the father swap, so need to set the age relative to its future starting point
		fixToTimes(oldfather,newdist);// sets the offset accounting for "inserting" a node

		// Update the fathers
		swapFather(notnewfather,revnode);// now the father of "notnewfather" is correct, but the father of which is which!
		swapFather(revnode,oldfather);//now all the fathers are correct

		swapEdgeTo(revnode,oldfather);
		revnode=getOldestReversedNode();
	}
	*whichto=orderNodes(which,dist,oldwhich); // reorder the nodes themselves
	updateList(which,*whichto,affedges);
	int numedgechanges=moveClonalFixFrom();//and update the "from" branches
	computeTTotal();
#if defined DEBUG
	try{testNodeAges(); //test its ok so far
	}catch(char * x) {
		dlog(1)<<x<<endl<<"Failed to move node "<<which<<" to "<<*whichto<<endl;
		for(int i=0;i<2*getN()-1;i++) { dlog(1)<<"Age["<<i<<"]="<<getNode(i)->getAge()<<endl;
	};exit(1);}
	testTree(); // this should be a valid tree at this stage
#endif
	return numedgechanges;
}

void RecTree::updateEdgeTimes(int e,double reldist)
{
	double newto=edge[e]->getTimeTo()*reldist;
	double newfrom=edge[e]->getTimeFrom()+newto-edge[e]->getTimeTo();
	edge[e]->setTimeTo(newto);
	edge[e]->setTimeFrom(newfrom);
}

int RecTree::moveClonalFixFrom(vector<int> ornodes)
{
    // We here decide what to do with an edge that is involved in the node move
    // This version keeps edges *relative* "to" position along the node constant
    int numedgechanges=0;
    bool useedges=false;
    if(ornodes.size()>0)
    {//use the vector of original nodes
        if((int)ornodes.size()!=numRecEdge())
        {
            cerr<<"moveClonalFixFrom: Wrong number of replacement edges!"<<endl;
            throw("moveClonalFixFrom: Wrong number of replacement edges!");
        }
        useedges=true;
    }
    for(unsigned int i=0;i<edge.size();i++)
    {
        // check that the "from" edge is still on the correct node
        while(getEdgeTimeAbsFrom(i) < getNode(edge[i]->getEdgeFrom())->getAge())
        {
            // choose a random daughter
            if(useedges)
            {// recall the daughter from ornodes
	        orderEdges(i);// give this edge the correct index *in the 0..i range*
                changeEdgeNodeFrom(i,ornodes[i]);// this has the wrong index.  we do this just to get a correctly sorted list at this stage
#if defined DEBUG
                if(getEdgeTimeAbsFrom(i) < getNode(edge[i]->getEdgeFrom())->getAge())
                {
		    cerr<<std::setprecision(64)<<newick(64)<<endl;
                    cerr<<std::setprecision(64)<<"Error in moveClonalFixFrom: edge "<<i<<" is supposed to be from node "<<ornodes[i]<<" but edges age "<<getEdgeTimeAbsFrom(i) <<" is < nodes age "<< getNode(edge[i]->getEdgeFrom())->getAge()<<endl;
                    throw("Edge Replacement error");
                }
#endif
                break;
            }
            else
            {
                if(gsl_rng_uniform(rng)<0.5)
                {
                    changeEdgeNodeFrom(i,getNode(edge[i]->getEdgeFrom())->getLeft()->getId());
                }
                else
                {
                    changeEdgeNodeFrom(i,getNode(edge[i]->getEdgeFrom())->getRight()->getId());
                }
                numedgechanges++;
            }
        }
        if(getNode(edge[i]->getEdgeFrom())->getId()!=root->getId())
        {
            while(getEdgeTimeAbsFrom(i) > getNode(edge[i]->getEdgeFrom())->getFather()->getAge())
            {
                changeEdgeNodeFrom(i,getNode(edge[i]->getEdgeFrom())->getFather()->getId());
                numedgechanges--;
                if(getNode(edge[i]->getEdgeFrom())->getFather()==NULL) break;
            }
        }
        orderEdges(i);// give this edge the correct index *in the 0..i range*
    }
    if(useedges)
    {// recall the daughter from ornodes, now in correct (from time) order
        for(unsigned int i=0;i<edge.size();i++)
        {
            changeEdgeNodeFrom(i,ornodes[i]);
        }
        return(-1);//don't return the number of changes if it is wrong!
    }
    return(numedgechanges);// return the net number of edge changes
}

void RecTree::orderEdges(int which)
{
    // keep swapping which with its neighbour until the ordering *UNTIL NODE which* is correct
    if(which>0)
    {
        while(getEdgeTimeAbsFrom(which-1)>getEdgeTimeAbsFrom(which))
        {
            swapEdge(which-1,which);
            which--;
            if(which==0)
                break;
        }
    }
}

void RecTree::swapEdge(int a, int b)
{
    RecEdge * efrom=edge[a];
    edge[a]=edge[b];
    edge[b]=efrom;
}

double RecTree::getEdgeTreeTime(int i) const
{
    int n1=edge[i]->getEdgeTo(),n2=edge[i]->getEdgeFrom();
    if(n1==n2)
        return(edge[i]->getTimeFrom()-edge[i]->getTimeTo());
    n1=nodes[n1]->getFather()->getId();
    double t1=nodes[n1]->getAge()-edge[i]->getTimeTo();
    double t2=0;
    if(n2==2*n-2)
    {
        t2=edge[i]->getTimeFrom();
    }
    else
    {
        n2=nodes[n2]->getFather()->getId();
        t2=nodes[n2]->getAge()-edge[i]->getTimeFrom();
    }
    while(n1!=n2)
    {
        if(nodes[n1]->getAge()<nodes[n2]->getAge())
        {
            n1=nodes[n1]->getFather()->getId();
            t1+=nodes[n1]->getDist();
        }
        else
        {
            n2=nodes[n2]->getFather()->getId();
            t2+=nodes[n2]->getDist();
        }
    }
    return(t1+t2);
}

void RecTree::updateList(int which, int newloc,vector<int> * list)
{
	for(unsigned int i=0;i<list->size();i++){
		if(list->at(i)==which) list->at(i)=newloc;
		else if(list->at(i)>which && list->at(i)<=newloc) list->at(i)--;
		else if(list->at(i)<which && list->at(i)>=newloc) list->at(i)++;
	}
}

int RecTree::orderNodes(int which, double dist,int oldwhich)
{
    int newloc=-1;
    // keep swapping which with its neighbour until the ordering is correct
    if(dist<0)
    {
        for(int i=which;i>=getN();i--)
        {
	    if(i==oldwhich||(getNode(i-1)->getAge()<getNode(i)->getAge() && (oldwhich<0||i<oldwhich))) {
                newloc=i;
                break;
	    }else{
                swapNode(i-1,i);
            }
        }
    }
    else
    {
        for(int i=which+1;i<2*getN()-1;i++)
        {
	    if(i-1==oldwhich||(getNode(i-1)->getAge()<getNode(i)->getAge() && (i-1>oldwhich ||oldwhich<0))) {
                newloc=i-1;
                break;
	    }else{
                swapNode(i-1,i);
            }
        }
    }
    if(newloc==-1)
        newloc=2*getN()-2; // if this is the oldest node it should stay where it is
    return newloc;
}

void RecTree::swapNode(int a, int b)
{
    Tree::swapNode(a,b);
    for(unsigned int e=0;e<edge.size();e++) {
	if(edge[e]->getEdgeTo()==a)       edge[e]->setEdgeTo(b);
	else if(edge[e]->getEdgeTo()==b)  edge[e]->setEdgeTo(a);
	if(edge[e]->getEdgeFrom()==a)     edge[e]->setEdgeFrom(b);
	else if(edge[e]->getEdgeFrom()==b)edge[e]->setEdgeFrom(a);
    }
}

void RecTree::moveEdges(int e1,int e2,double t1,double t2)
{
	if(t2<0.0 && e1!=root->getId()) t2=getNode(e1)->getFather()->getAge();
	for(unsigned int i=0;i<edge.size();i++){
		if(edge[i]->getEdgeTo()==e1 && getEdgeTimeAbsTo(i)>t1 &&getEdgeTimeAbsTo(i)<t2) edge[i]->setTimeTo(getEdgeTimeAbsTo(i)-getNode(e2)->getAge(),e2);
		if(edge[i]->getEdgeFrom()==e1 && getEdgeTimeAbsFrom(i)>t1 &&getEdgeTimeAbsFrom(i)<t2) edge[i]->setTimeFrom(getEdgeTimeAbsFrom(i)-getNode(e2)->getAge(),e2);
	}
}

vector<int> RecTree::getAffEdges(vector<int> * samplespace)
{
	vector<int> validedges;
	int eto,efrom;
	for(int e=0;e<numRecEdge();e++){
		eto = getEdge(e)->getEdgeTo();
		efrom = getEdge(e)->getEdgeFrom();
		for(unsigned int i=0;i<samplespace->size();i++) {
			if(eto==samplespace->at(i)||efrom==samplespace->at(i)) {
				validedges.push_back(e);
				i=samplespace->size();
			}
		}
	}
	return(validedges);
}

int RecTree::sampleEdge(vector<int> * samplespace)
{
	if(numRecEdge()==0) return(-1);// no valid edges
	int which=gsl_rng_uniform_int(rng,numRecEdge());
	if(which>=numRecEdge()) throw("Selected invalid edge!");
	if(samplespace==NULL) return(which);
	vector<int> validedges = getAffEdges(samplespace);
	if(validedges.size()==0) return(-1);// no valid edges
	int whichve=gsl_rng_uniform_int(rng,validedges.size());
	if(whichve>=(int)validedges.size()) throw("Selected edge not in samplespace!");
	which=validedges[whichve];
	if(which>=numRecEdge()) throw("Selected invalid edge!");
	return(which);
}

int RecTree::lastCommonAncestor(int s1,int s2)
{
	// assume that the local tree is correct
	vector<int> r1,r2;
	int lcaindexr1=-1;/// index of last common ancestor in list from s1
	r1.push_back(tabFather[s1]);
	while(r1.back()>=0) {
		r1.push_back(tabFather[r1.back()]);
	}
	r2.push_back(tabFather[s2]);
	while(r2.back()>=0) {
		r2.push_back(tabFather[r2.back()]);
	}
	for(unsigned int i=1;i<=max(r1.size(),r2.size());i++) {
		if(r1[r1.size()-i]==r2[r2.size()-i]) {
			lcaindexr1=r1.size()-i;
		}else break;
	}
	return(r1[lcaindexr1]);
}

double RecTree::pairwiseDist(int s1,int s2)
{
	int lca=lastCommonAncestor(s1,s2);
	return(2.0*age[lca]);
}

vector<vector<double> > RecTree::pairwiseDistanceMatrix(long site)
{
	vector<vector<double> > res;
	makeLocalTreeKeepingFathers(site);
	for(int i=0;i<getN();i++) res.push_back(vector<double>(getN(),0.0));
	for(int i=0;i<getN();i++) {
		for(int j=i+1;j<getN();j++) {
			res[i][j]=res[j][i]=pairwiseDist(i,j);
		}
	}
	return(res);
}




} // end namespace weakarg
