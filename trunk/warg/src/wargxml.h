#ifndef WARGXML_H
#define WARGXML_H
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ios>

namespace weakarg
{

/**
    @brief Reads in the output file
*/
class WargXml
{
public:
	WargXml(std::string f);
	~WargXml();
	std::streampos getLineContaining(std::string str,bool getlast=false);
	std::streampos gotoLineContaining(std::string str,bool getlast=false);
	inline std::streampos getLastIteration(){
		return(getLineContaining("<Iteration>",true));
	};///<Gets the last iteration
	inline std::string getLine(){
		std::string res; 
		getline(iterfile,res); 
		return(res);
	}
	inline void seekg(std::streampos is){iterfile.seekg(is);};
	inline std::streampos tellg(){return(iterfile.tellg());};
	std::string getTree(std::streampos itstart=-1,bool isfinal=true);// if isfinal, gets the last iteration's tree; else gets the next
	std::string getParam(std::string tag, std::streampos itstart=-1);
	inline bool isempty(){
		bool ret=false;std::streampos sp=tellg(); restart();
		if(iterfile.peek()==EOF) ret=true;
		seekg(sp);
		return(ret);
	}
	inline bool eof(){return(iterfile.eof());};
	inline void restart(){iterfile.clear();iterfile.seekg(std::ios_base::beg);}
	inline void clear(){iterfile.clear();};
private:
	std::string fname;
	std::ifstream iterfile;
};

}//end namespace
#endif
