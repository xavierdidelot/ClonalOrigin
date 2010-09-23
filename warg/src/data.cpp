#include "data.h"
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;
namespace weakarg
{

Data::Data(int n,vector<int> blocks)
{
    this->n=n;
    this->L=blocks.back();
    this->blocks=blocks;
    data=vector<string>(n,string(L,'N'));
    poly=vector<bool>(L,true);
    begEnd=vector<bool>(L+1,false);
    for (unsigned int i=0;i<blocks.size();i++)
        begEnd[blocks[i]]=true;
   for (unsigned int i=0;i<blocks.size()-1;i++) regions.push_back(i);
 
}

Data::Data(string filename)
{
    string line;
    ifstream file;
    file.open(filename.data());//Open file
    unsigned int which=-1;
    while (1)
    {
        getline(file,line);//Read next line from file
        if (file.eof())
            break;//Stop if end of file
	if (line[0]=='=') which=-1;
        if (line.size()==0 || line[0]=='#' || line[0]=='=')
            continue;//Ignore empty lines, comments, and end of block lines
        if (line[0]=='>')
        {//Header line
            line.erase(0,1);
	    which++;
	    string tmpname=cleanName(line);
	// aught to check the name matches if block !=0, but we don't at the moment
            if (which>=data.size()){
		addIsolate(tmpname);
	    }else if(getNumber(tmpname)<0){
	            cerr<<"Data names are inconsistent: "<<tmpname<<" not found!"<<endl;
            throw("Data error");
	    }else if(getNumber(tmpname)!=which){
	            cerr<<"Data names are inconsistent: "<<tmpname<<" is not the "<<which<<"th name in the first block!"<<endl;
            throw("Data error");
	    }
            if (which==0)
            {
                blocks.push_back(data[0].size());
            };
            continue;
        }
        //Sequence data line
        data[which].append(line);
    }
    file.close();//Close file
    n=data.size();
    L=data[0].size();
    blocks.push_back(L);
    for (unsigned int i=1;i<n;i++)
        if (data[0].size()!=data[i].size())
        {
            cerr<<"Data is inconsistent: "<<data[0].size()<<"!="<<data[i].size()<<endl;
            break;
        }
    dlog(1)<<"Read input file with "<<n<<" isolates and "<<getB()<<" blocks for a total of "<<L<<" sites."<<endl;
    for (unsigned int i=0;i<n;i++)
        for (unsigned int j=0;j<L;j++)
            data[i][j]=convert(data[i][j]);
    poly=vector<bool>(L,true);
    for (unsigned int i=0;i<L;i++)
        makePoly(i);
    begEnd=vector<bool>(L+1,false);
    for (unsigned int i=0;i<blocks.size();i++)
        begEnd[blocks[i]]=true;
   for (unsigned int i=0;i<blocks.size()-1;i++) regions.push_back(i);
}

char Data::convert(char in)
{
    switch (in)
    {
    case 'a':
    case 'A':
        return 0;
    case 't':
    case 'T':
        return 1;
    case 'c':
    case 'C':
        return 2;
    case 'g':
    case 'G':
        return 3;
    default:
        return 'N';
    }
}

char Data::convertBack(char in)
{
    switch (in)
    {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    default:
        return 'N';
    }
}

void Data::output(ostream * out)
{
    for (int i=0;i<getB();i++)
    {
        for (unsigned int j=0;j<n;j++)
        {
            *out<<">"<<j<<endl;
            for (int k=blocks[i];k<blocks[i+1];k++)
                *out<<convertBack(data[j][k]);
            *out<<endl;
        }
        *out<<"="<<endl;
    }
}

string Data::cleanName(string str)
{
	const string match="\">;, ";
	size_t found=str.find_first_of(match);
	while(found!=std::string::npos) {
		str.erase(found,1);
		found=str.find_first_of(match);
	}
	return(str);	
}

double Data::watterson()
{
    int p=0;
    for (unsigned int i=0;i<poly.size();i++)
        if (poly[i])
            p++;//bad with missing data
    double s=0;
    for (unsigned int i=1;i<=n;i++)
        s+=1.0/i;
    return 1.0*p/s;
}

Data::~Data()
{}

void Data::subset(vector<int> numv,int seed)
{
  if(numv.size()==0) return;
  int num=numv[0];
  if (num<=0 &&seed>=0 || num<0) return;
  if (num>getB() &&seed>=0) {cerr<<"Not enough regions for subset, ignoring argument."<<endl;return;}
  if(seed>=0)dlog(1)<<"Picking subset of "<<num<<" regions using seed "<<seed<<"..."<<endl;
  else dlog(1)<<"Using specified subset of regions..."<<endl;
  //Choose regions to include
  regions.clear();
  vector<bool> tried(getB(),false);
  if(seed<0) {
    for(int i=0;i<numv.size();i++) if(numv[i]<tried.size()&&numv[i]>=0) {if(!tried[numv[i]]) {regions.push_back(numv[i]);tried[numv[i]]=true;}else{cerr<<"Warning: Duplicate blocks of "<<numv[i]<<" were specified and have been removed from the analysis."<<endl;}}else{cerr<<"Warning: Invalid blocks were specified and have been removed from the analysis."<<endl;}
    num=regions.size();
  }else{
  for (int i=0;i<num;i++) {
    seed=((seed+1)*99)%getB();
    while (tried[seed]) seed=(seed+1)%getB();
    tried[seed]=true;
    regions.push_back(seed);
  }
  }
  //Modify data accordingly
  for (int i=0;i<getN();i++) {
    string tmp;
    for (int j=0;j<num;j++) {
    string tmp2=data[i].substr(blocks[regions[j]],blocks[regions[j]+1]-blocks[regions[j]]);
    tmp.append(tmp2);
   }
  data[i]=tmp;
  }
  //Modify blocks accordingly
  vector<int> blocks2;
  blocks2.push_back(0);
  for (int j=0;j<num;j++) blocks2.push_back(blocks2.back()+blocks[regions[j]+1]-blocks[regions[j]]);
  blocks=blocks2;
  L=blocks.back();
  //Recompute poly
    poly=vector<bool>(L,true);
  for (unsigned int i=0;i<L;i++)
      makePoly(i);
  //Recompute bedEnd
  begEnd=vector<bool>(L+1,false);
  for (unsigned int i=0;i<blocks.size();i++)
      begEnd[blocks[i]]=true;
  dlog(1)<<"Data has now "<<getB()<<" blocks for a total of "<<L<<" sites."<<endl;
}


void Data::readRegionsFromFile(WargXml * infile)
{
	streampos sp=infile->tellg();
	infile->restart();
	string sregions=infile->getParam("regions",-1);
	vector<int> iregions;
	size_t pos = 0, ppos = -1;
	    while (pos != string::npos) {
		pos=sregions.find_first_of(',',ppos+1);
		int reg=atoi(sregions.substr(ppos+1,pos-ppos-1).c_str());
		iregions.push_back(reg);
		ppos=pos;
  	    }
	subset(iregions,-1);
	infile->seekg(sp);
}


} // end namespace weakarg
