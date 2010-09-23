#include "data.h"
#include <fstream>
#include <sstream>

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
}

Data::Data(string filename)
{
    string line;
    ifstream file;
    file.open(filename.data());//Open file
<<<<<<< .mine
    if(!file.is_open())
    {
    	cerr << "Unable to read file \""<<filename<<"\".\n";
    	throw "File unreadable\n";
    }
    unsigned int which=0;
=======
    unsigned int which=-1;
>>>>>>> .r336
    blocks.push_back(0);
    while (1)
    {
        getline(file,line);//Read next line from file
        if (file.eof())
            break;//Stop if end of file
<<<<<<< .mine
        if(line[0]=='=')
        {
        	// pad all sequences with - to match current size
        	// this enables the subset blocks in Mauve XMFA format to be read
        	unsigned int maxblock = 0;
        	for(unsigned int i=0; i<data.size(); i++){
        		if(maxblock < data[i].size()) maxblock = data[i].size();
        	}
		data.resize(29);
        	for(unsigned int i=0; i<29; i++){
        		data[i].resize(maxblock, '-');
        	}
            blocks.push_back(data[0].size());
        }
=======
	if (line[0]=='=') which=-1;
>>>>>>> .r336
        if (line.size()==0 || line[0]=='#' || line[0]=='=')
            continue;//Ignore empty lines, comments, and end of block lines
        if (line[0]=='>')
        {//Header line
            line.erase(0,1);
<<<<<<< .mine
            istringstream iss(line);
            iss>>which;
            which--;	// numbering starts at 1 in XMFA, indices here start at 0
            while (which>=data.size())
                data.push_back("");
=======
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
>>>>>>> .r336
            continue;
        }
        //Sequence data line
        data[which].append(line);
    }
    file.close();//Close file
    n=data.size();
    L=data[0].size();
//    blocks.push_back(L);
    for (unsigned int i=1;i<n;i++)
        if (data[0].size()!=data[i].size())
        {
            cerr<<"Data is inconsistent: "<<data[0].size()<<"!="<<data[i].size()<<endl;
            break;
        }
    cout<<"Read input file with "<<n<<" isolates and "<<getB()<<" blocks for a total of "<<L<<" sites."<<endl;
    for (unsigned int i=0;i<n;i++)
        for (unsigned int j=0;j<L;j++)
            data[i][j]=convert(data[i][j]);
    poly=vector<bool>(L,true);
    for (unsigned int i=0;i<L;i++)
        makePoly(i);
    begEnd=vector<bool>(L+1,false);
    for (unsigned int i=0;i<blocks.size();i++)
        begEnd[blocks[i]]=true;
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


} // end namespace weakarg
