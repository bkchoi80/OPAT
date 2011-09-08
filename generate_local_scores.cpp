#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <stdio.h>

#include "opat.hpp"
#include "generate_local_scores.hpp"


using namespace std;

void generateLocalScores(const vector< vector<double> > &data, size_t maxInDegree, size_t numLookAheads, double rho)
{
	vector<size_t> family;
	list<size_t> candidates;
	vector<ofstream *> files;
	for(size_t nodeNum = 0; nodeNum < data[0].size(); nodeNum++) {
		candidates.push_back(nodeNum);
		ostringstream filename;
		filename << "scores_node" << nodeNum << ".txt";
		files.push_back(new ofstream(filename.str().c_str()));
	}

	for(size_t familySize = 1; familySize <= min(maxInDegree + 1, candidates.size()); familySize++)
		generateLocalScoresSub(files, familySize, candidates.begin(), candidates.size(), family, data, numLookAheads, rho);

	for(size_t nodeNum = 0; nodeNum < data[0].size(); nodeNum++)
		files[nodeNum]->close();
}

void generateLocalScoresSub(vector<ofstream *> &files, size_t numAdditions, list<size_t>::const_iterator candidates, size_t numCandidates, vector<size_t> &family, const vector< vector<double> > &data, size_t numLookAheads, double rho)
{
	if(numCandidates < numAdditions)
		return;

	if(numAdditions == 0) {
		OPAT tree(data, family);
		tree.build(numLookAheads, rho);
		for(size_t familyNum = 0; familyNum < family.size(); familyNum++)
			(*files[family[familyNum]]) << tree.computeConditionalML(familyNum) << endl;
	}
	else {
		family.push_back(*candidates);
		candidates++;
		generateLocalScoresSub(files, numAdditions - 1, candidates, numCandidates - 1, family, data, numLookAheads, rho);
		family.pop_back();
		generateLocalScoresSub(files, numAdditions, candidates, numCandidates - 1, family, data, numLookAheads, rho);
		candidates--;
	}
}
