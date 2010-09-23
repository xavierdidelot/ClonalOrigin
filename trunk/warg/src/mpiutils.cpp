/*
 * mpiutils.cpp
 *
 *  Created on: Mar 2, 2009
 *      Author: koadman
 */

#include "mpiutils.h"
#include "weakarg.h"
#include <sstream>
#include <fstream>

using namespace std;

namespace weakarg {

int mpirank;
int mpinodecount;
std::string logfilename;

#ifdef WEAKARG_MPI

void initmpi( int argc, char** argv )
{
	// Set up the MPI environment
    MPI_Init(&argc, &argv);	// some MPI options can be passed on the command-line
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);	// get these values from the mpi library
    MPI_Comm_size(MPI_COMM_WORLD, &mpinodecount);
}

void endmpi()
{
	MPI_Finalize();
}

ostream& dlog()
{
	static ostream* out = NULL;
	if(out==NULL)
	{
		if(logfilename.length() == 0)
		{
			out = &cout;
		}else{
			stringstream pname;
			pname << logfilename << "." << mpirank;
			out = new ofstream();
			dynamic_cast<ofstream*>(out)->open( pname.str().c_str() );
			if(!dynamic_cast<ofstream*>(out)->is_open())	cerr << mpirank << " error opening log file " << pname.str() << endl;
		}
	}
	return *out;
}

mpiofstream::mpiofstream(const std::string& fname) :
	fname(fname)
{
	int rval = MPI_File_open(MPI_COMM_WORLD, (char*)fname.c_str(), MPI_MODE_WRONLY | MPI_MODE_UNIQUE_OPEN | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
	cerr << "Opening file " << fname << " returned " << rval << endl;
}
void mpiofstream::flush()
{
	MPI_Status status;
	char* buf;
	MPI_Alloc_mem(ss.str().length()+1, MPI_INFO_NULL, &buf);
	strcpy(buf, ss.str().c_str());
	MPI_File_write_shared(fh, buf, ss.str().length(), MPI_CHAR, &status);
	MPI_Free_mem(buf);
	ss.str("");
}
mpiofstream::~mpiofstream()
{
	MPI_File_close(&fh);
}
std::ostream& mpiofstream::operator()(bool master)
{
	static nullstream nullstr;
	if(master && mpirank!=0)	return nullstr;
	return ss;
}


#else	// serial implementations

void initmpi( int argc, char** argv ){}
void endmpi(){}
ostream& dlog()
{
	static ostream* out = NULL;
	if(out==NULL)
	{
		if(logfilename.length() == 0)
		{
			out = &cout;
		}else{
			stringstream pname;
			pname << logfilename;
			out = new ofstream();
			dynamic_cast<ofstream*>(out)->open( pname.str().c_str() );
			if(!dynamic_cast<ofstream*>(out)->is_open())	cerr << " error opening log file " << pname.str() << endl;
		}
	}
	return *out;
}

mpiofstream::mpiofstream(const std::string& fname) :
	fname(fname)
{
	file = new ofstream(fname.c_str());
}
void mpiofstream::flush()
{
	*file << ss.str();
	ss.str("");
}
mpiofstream::~mpiofstream(){ delete file; }
std::ostream& mpiofstream::operator()(bool master)
{
	return ss;
}

#endif	// WEAKARG_MPI


//
// Functions below are identical regardless of MPI
//


ostream& dlog(int verbosity)
{
	static nullstream nullstr;
	if(verbosity > opt().verbose )
	{
		return nullstr;
	}else{
		return dlog();
	}
}



} // end namespace weakarg
