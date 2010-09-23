/*
 * weakarg.h
 *
 *  Created on: Mar 2, 2009
 *      Author: koadman
 */

#ifndef WEAKARG_H_
#define WEAKARG_H_

#include <vector>
#include <string>

namespace weakarg
{

static const int NUMMOVES=11;

class ProgramOptions
{
public:
	ProgramOptions() :
		theta(-1.0),
		rho(-1.0),
		thetaPerSite(false),
		rhoPerSite(false),
		delta(-1.0),
		preburnin(100000),
		burnin(100000),
		additional(100000),
		thinin(100),
		allowclonal(true),
		temperreps(0),
		temperT(1.0),
		subset(std::vector<int>(0)),
		subsetSeed(-1),
		greedyWeight(0),
		verbose(0)
		{
			movep.resize(NUMMOVES,1.0);
		}
    std::vector<double> movep; //length (NUMMOVES) vector with weights set to default 1

	double theta;
	double rho;
	bool thetaPerSite;
	bool rhoPerSite;
	double delta;
	int preburnin;
	int burnin;
	int additional;
	int thinin;
	bool allowclonal;
	int temperreps;
	double temperT;
	std::vector<int> subset;
	int subsetSeed;
	double greedyWeight;
	int verbose;
	std::string outfile;
	std::string csvoutfile;
	std::string logfile;
};
ProgramOptions& opt();

} // end namespace weakarg

#endif /* WEAKARG_H_ */
