#include <string>
#include <iostream>
#include <stdlib.h>

#include "myutil.hpp"
#include "generate_local_scores.hpp"


using namespace std;

int main(int argc, char *argv[]){
	string datafile;
	vector< vector<double> > data;
	if(argc < 2){
    	cout << "Data file: ";
    	cin >> datafile;
    }
    else
    	datafile = string(argv[1]);
	if(!myutil::readDelim(datafile.c_str(), data))
		exit(-1);

	size_t maxInDegree;
	size_t numLookAheads;
	string outputfile;
    cout << "Maximum in-degree: ";
    cin >> maxInDegree;
    cout << "Number of lookaheads: ";
    cin >> numLookAheads;

    generateLocalScores(data, maxInDegree, numLookAheads);
}

