/*
 * mpiutils.h
 * MPI implementations of functions go here, along with serial implementations
 * Hopefully we can localize the #ifdef #else #endif statements to this file
 *
 *  Created on: Mar 2, 2009
 *      Author: koadman
 */

#ifndef MPIUTILS_H
#define MPIUTILS_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#ifdef WEAKARG_MPI
#include <mpi.h>
#endif	// WEAKARG_MPI

namespace weakarg
{

extern int mpirank;	///< The rank (process ID) of the running mpi thread
extern int mpinodecount;	///< The total number of MPI processes
extern std::string logfilename;

void initmpi( int argc, char** argv );	///< Set up the MPI environment
void endmpi();	///< Tear down the MPI environment

std::ostream& dlog();	///< return a logging output stream

std::ostream& dlog(int verbosity);	///< return a logging output stream dependent on a verbosity level

/*** A class that sends output directly to the bitbucket */
struct nullstream:
	std::ostream {
	nullstream(): std::ios(0), std::ostream(0) {}
};

/**
 * A class to write either to a standard file or a shared file via MPI-IO.
 * After performing one or more writes with the << operator, flush() must be called
 * to force the write.  This class merely buffers data in memory until flush() is called.
 * In the non-MPI case, this class will forward all writes to a standard file.
 * In the MPI case, the class will use shared I/O to append data to the file.
 */
class mpiofstream
{
public:
	mpiofstream(const std::string& fname);///<Opens a file with the given name
	~mpiofstream();
	void flush();///<Causes data to be written
	std::ostream& operator()(bool master);///<returns an output stream.  if master is true, then only the master MPI process will be allowed to write.

protected:
	std::string fname;
	std::ostringstream ss;
#ifdef WEAKARG_MPI
	MPI_File fh;
#else
	std::ofstream* file;
#endif
};


}	// end namespace weakarg

#endif /* MPIUTILS_H_ */
