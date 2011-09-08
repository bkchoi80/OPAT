#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

#include "myutil.hpp"


namespace myutil {
	using namespace std;

	bool readDelim(const char* const filename, vector< vector<double> >& data)
	{
		ifstream ifs(filename);
		if(!ifs) {
			cerr << "Cannot read the file. Parse failed." << endl;
			return(false);
		}
		data.clear();

		stringstream iss;
		string line;
		double element;

		while(getline(ifs, line)) {
			if(line.size() == 0)
				break;

			data.resize(data.size() + 1);
			iss << line;
			while(iss >> element)
				data.back().push_back(element);

			if((data.size() > 1) & (data.back().size() != data[data.size() - 2].size())) {
				cerr << "The size of line " << data.size() - 1 << " and line " << data.size() << " are different. Parse failed." << endl;
				return(false);
			}
			iss.clear();
		}

		return(true);
	}
}
