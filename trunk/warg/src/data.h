#ifndef DATA_H
#define DATA_H
#include <vector>
#include <string>
#include <iostream>
#include "mpiutils.h"
#include "wargxml.h"

namespace weakarg
{

/**
    @brief Sequence data for the leaves of the topology
*/
class Data
{
public:
    Data(std::string filename);///<Reads in the data from a file
    Data(int n,std::vector<int> blocks);///<Creates empty data
    ~Data();
    inline char get(int i,int j)
    {
        return data[i][j];
    } ///<Get accessor to the data
    inline void set(int i,int j,char c)
    {
        data[i][j]=c;
        makePoly(j);
    } ///<Set accessor to the data
    inline void set_NO_POLY_UPDATE(int i,int j,char c)
    {
        data[i][j]=c;
    }
    inline int getN()
    {
        return n;
    }///<Returns the number of isolates
    inline int getL()
    {
        return L;
    }///<Returns the length of the sequences
    inline int getB()
    {
        return blocks.size()-1;
    }///<Returns the number of blocks
    inline int inblock(int site)
    {
        int blockin=0;
        while(blocks[blockin+1]<=site)
            blockin++;
        return blockin;
    }///<Returns the block a given site is in
    inline std::vector<int> * getBlocks()
    {
        return &blocks;
    }///<Returns the block structure
    void output(std::ostream * out);
    inline bool isPoly(unsigned int site)
    {
        return poly[site];
    }///<Returns whether a site is polymorphic
    void makePoly(unsigned int site)
    {
        if (data[0][site]>3)
        {
            poly[site]=true;
            return;
        }
        for (unsigned int i=1;i<n;i++)
            if (data[0][site]!=data[i][site])
            {
                poly[site]=true;
                return;
            }
        poly[site]=false;
    }
    double watterson();///<Returns watterson's estimate of theta
    inline int numPoly()
    {
        int r=0;
        for (unsigned int i=0;i<L;i++)
            if (poly[i])
                r++;
        return r;
    }///<Returns the number of polymorphic sites
    inline bool isBegEnd(int a)
    {
        return begEnd[a];
    }
    std::string cleanName(std::string tname);///< Cleans a name in a standard way
    inline void printNames(std::ostream& out){
	for(unsigned int i=0;i<datanames.size();i++) {
	out<<i<<","<<datanames[i]<<";";
	}
    }///< Prints an index->name mapping
    inline int getNumber(std::string str){
	str=getSname(cleanName(str));
	for(unsigned int c1=0;c1<datasnames.size();c1++){
		if(datasnames[c1].compare(str)==0) return(c1);
	}
	return(-1);
    }///<Gets the index of the isolate name specified (-1 if not found)
    inline void addIsolate(std::string str){
	datanames.push_back(str);
	datasnames.push_back(getSname(str));
	data.push_back("");
    }///<Adds an islate including its data, name and number 
    inline std::string getSname(std::string str){
	size_t found=str.find_first_of(":±");
	if(found!=std::string::npos) str.erase(found);
	return(str);
    }
    void subset(std::vector<int> numv,int seed);///<Use only a subset of the regions of size "num" and determined by seed "seed"
    std::vector<int> * getRegions(){return(&regions);}///<Gets the regions used in the analysis
    void readRegionsFromFile(WargXml * infile);///<reads the regions in from a file
protected:
    char convert(char in);///<Converts a character from A,C,G,T to 0,1,2,3
    char convertBack(char in);///<Converts a character from 0,1,2,3 to A,C,G,T
    unsigned int n;///<Number of isolates
    unsigned int L;///<Length of concatenated sequences
    
    std::vector<std::string> data;///<Concatenated data
    std::vector<std::string> datanames;///<Names of the data
    std::vector<std::string> datasnames;///<Short names of the data: numbers if using correct xmfa
    std::vector<int> blocks;///<Starting points of blocks in concatenated data, finished with L
    std::vector<int> regions;///<List of the regions used in this run
    
    std::vector<bool> poly;///<Indicates whether the sites are polymorphic
    std::vector<bool> begEnd;///<Indicates whether the sites are beginning or end of blocks
};

} // end namespace weakarg
#endif
