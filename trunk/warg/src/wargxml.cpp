#include "wargxml.h"
#include <fstream>
#include <sstream>
#include <cstdlib>

using namespace std;
namespace weakarg
{

WargXml::WargXml(std::string f)
{
	fname=f;
	iterfile.open(fname.data());
        if (!iterfile) {
            cerr << "Can't open input file " << f << endl;
            exit(1);
        }
}

WargXml::~WargXml()
{
	if(iterfile.is_open()) iterfile.close();
}

streampos WargXml::gotoLineContaining(std::string str,bool getlast){
	streampos pos =iterfile.tellg (),it=-1;
	string res;
	size_t found;
	int len=str.length();
	while(1){
		pos =iterfile.tellg ();
		getline(iterfile,res);
		found=res.find_first_of('<');
		if(found!=string::npos) if(res.length()>=len+found) if(res.substr(found,len).compare(str)==0) {it=pos;}
		if (iterfile.eof()||!iterfile.good()||(!getlast && it>=0)) break;
	}
	iterfile.clear();
	iterfile.seekg(it);
	return(it);
}

std::streampos WargXml::getLineContaining(std::string str,bool getlast){
	std::streampos p0=iterfile.tellg();
	std::streampos r=gotoLineContaining(str,getlast);
	iterfile.seekg(p0);
	return(r);
}

std::string WargXml::getParam(string tag, streampos itstart)
{
	string res;
	string stag=tag,etag=tag;
	stag.insert(0,"<");
	etag.insert(0,"</");
	stag.append(">");
	etag.append(">");
	if(itstart>=0) iterfile.seekg(itstart);
	else itstart=iterfile.tellg ();
	streampos p=gotoLineContaining(stag.c_str(),false);
	if(p<0){
		cout<<"Tag "<<stag<<" not found in "<<fname<<endl;
		throw("Tag not found in file!");
	}
	getline(iterfile,res);
	size_t found=res.find(stag.c_str());
	res.erase(found,stag.length());
	size_t found2=res.find(etag.c_str());
	res=res.substr(found,found2);
	iterfile.seekg(itstart);
	return(res);
}

std::string WargXml::getTree(streampos itstart,bool isfinal)
{
	string res;
	size_t found;
	if(itstart>=0) iterfile.seekg(itstart);
	else itstart=iterfile.tellg ();
	getLineContaining("<Tree>",isfinal);
	while((found=res.find_first_of('('))==string::npos) {
		getline(iterfile,res);
		if(iterfile.eof())throw("Error in readtree: <Tree> does not appear to contain a newick tree!");
	}
	size_t found2=res.find("</Tree>");
	if(found2!=string::npos) found2-=6;
	res=res.substr(found,found2);
	iterfile.seekg(itstart);
	return(res);
}

}//end namespace
