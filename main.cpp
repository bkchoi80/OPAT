#include <vector>
#include <math.h>

#include "opat.hpp"
#include "myutil.hpp"
#include "multinomial.hpp"

using namespace std;

int main(int argc, char *argv[]) {
	string datafile;

	if (argc < 2) {
		cout << "Data file: ";
		cin >> datafile;
	} else
		datafile = string(argv[1]);

	vector< vector<double> > data;
	myutil::readDelim(datafile.c_str(), data);

	vector<size_t> type(data[0].size());
	for(size_t varNum = 0; varNum < data[0].size(); varNum++)
		type[varNum] = floor(data[0][varNum] + 0.5);
	data.erase(data.begin());



	// cout << endl << "Initializing...";
	// OPAT tree(data);

	double rho;
	cout << endl << endl << "Rho: ";
	cin >> rho;

	size_t numLookAhead;
	cout << "Number of Lookahead: ";
	cin >> numLookAhead;

	// cout << endl << "Building tree...";
	// tree.build(numLookAhead, rho);
	// cout << endl;

	// tree.debug();
	// cout << endl;





	std::ofstream ofs;
	ofs.open("output.txt", std::ios_base::trunc);

	vector< vector<double> > MI(data[0].size());
	for(size_t dimNum = 0; dimNum < data[0].size(); dimNum++)
		MI[dimNum] = vector<double> (data[0].size());
	vector< vector<double> > subdata1(data.size());
	for(size_t dataNum = 0; dataNum < subdata1.size(); dataNum++)
		subdata1[dataNum] = vector<double>(1);
	vector< vector<double> > subdata2(data.size());
	for(size_t dataNum = 0; dataNum < subdata2.size(); dataNum++)
		subdata2[dataNum] = vector<double>(2);
	vector<double> partialML(data[0].size());

	for(size_t dimNum = 0; dimNum < data[0].size(); dimNum++) {
		for(size_t dataNum = 0; dataNum < subdata1.size(); dataNum++)
			subdata1[dataNum][0] = data[dataNum][dimNum];
		if(type[dimNum] == 0) {
			OPAT tree(subdata1);
			tree.build(numLookAhead, rho);
			partialML[dimNum] = tree.computeML();
		}
		else {
			partialML[dimNum] = multinomial::computeLogLikelihood(subdata1);
		}
	}

	for(size_t dimNum1 = 0; dimNum1 < data[0].size(); dimNum1++) {
		for(size_t dimNum2 = 0; dimNum2 < dimNum1; dimNum2++) {
			if(type[dimNum1] == 0) {
				if(type[dimNum2] == 0) {
					for(size_t dataNum = 0; dataNum < subdata2.size(); dataNum++) {
						subdata2[dataNum][0] = data[dataNum][dimNum1];
						subdata2[dataNum][1] = data[dataNum][dimNum2];
					}
					OPAT tree(subdata2);
					tree.build(numLookAhead, rho);
					MI[dimNum1][dimNum2] = tree.computeML();
				}
				else {
					for(size_t valueNum = 0; valueNum < type[dimNum2]; valueNum++) {
						vector< vector<double> > filteredData;
						for(size_t dataNum = 0; dataNum < data.size(); dataNum++)
							if(data[dataNum][dimNum2] == valueNum)
								filteredData.push_back(vector<double>(1, data[dataNum][dimNum1]));
						OPAT tree(filteredData);
						tree.build(numLookAhead, rho);
						if(filteredData.size() > 0)
							MI[dimNum1][dimNum2] += tree.computeML();
					}
					MI[dimNum1][dimNum2] += partialML[dimNum2];
				}
			}
			else {
				if(type[dimNum2] == 0) {
					for(size_t valueNum = 0; valueNum < type[dimNum1]; valueNum++) {
						vector< vector<double> > filteredData;
						for(size_t dataNum = 0; dataNum < data.size(); dataNum++)
							if(data[dataNum][dimNum1] == valueNum)
								filteredData.push_back(vector<double>(1, data[dataNum][dimNum2]));
						OPAT tree(filteredData);
						tree.build(numLookAhead, rho);
						if(filteredData.size() > 0)
							MI[dimNum1][dimNum2] += tree.computeML();
					}
					MI[dimNum1][dimNum2] += partialML[dimNum1];
				}
				else {
					for(size_t dataNum = 0; dataNum < subdata2.size(); dataNum++) {
						subdata2[dataNum][0] = data[dataNum][dimNum1];
						subdata2[dataNum][1] = data[dataNum][dimNum2];
					}
						MI[dimNum1][dimNum2] = multinomial::computeLogLikelihood(subdata2);
				}
			}

			MI[dimNum1][dimNum2] -= partialML[dimNum1] + partialML[dimNum2];
			MI[dimNum1][dimNum2] /= data.size();
		}
	}

	for(size_t dimNum1 = 0; dimNum1 < data[0].size(); dimNum1++) {
		ofs << MI[dimNum1][0];
		for(size_t dimNum2 = 1; dimNum2 < data[0].size(); dimNum2++)
				ofs << "\t" << MI[dimNum1][dimNum2];
		ofs << endl;
	}

	ofs.close();
}

