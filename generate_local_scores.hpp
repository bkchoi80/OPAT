#ifndef GENERATE_LOCAL_SCORES_HPP_
#define GENERATE_LOCAL_SCORES_HPP_

#include <vector>
#include <list>

void generateLocalScores(const std::vector< std::vector<double> > &data, size_t maxInDegree, size_t numLookAheads = 2, double rho = 0.8);
void generateLocalScoresSub(std::vector<std::ofstream *> &files, size_t numAdditions, std::list<size_t>::const_iterator candidates, size_t numCandidates, std::vector<size_t> &family, const std::vector< std::vector<double> > &data, size_t numLookAheads, double rho);

#endif /* GENERATE_LOCAL_SCORES_HPP_ */
