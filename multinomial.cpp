/*
 * multinomial.cpp
 *
 *  Created on: May 25, 2010
 *      Author: bkchoi
 */

#include "multinomial.hpp"
#include <math.h>
#include <iostream>


using namespace std;

double multinomial::computeLogLikelihood(const vector< vector<double> > &data)
{
	double ll = 0;

	if(data[0].size() == 1) {
		vector<double> countTable;
		for(size_t dataNum = 0; dataNum < data.size(); dataNum++) {
			size_t newDatum = floor(data[dataNum][0] + 0.5);
			if(newDatum >= countTable.size())
				countTable.resize(newDatum + 1);
			countTable[newDatum]++;
		}

		for(size_t cellNum = 0; cellNum < countTable.size(); cellNum++)
			ll += countTable[cellNum] * (log(countTable[cellNum] + 1.0) - log((double) data.size() + (double) countTable.size()));
	}
	else {
		vector< vector<double> > countTable;
		for(size_t dataNum = 0; dataNum < data.size(); dataNum++) {
			size_t newDatum0 = floor(data[dataNum][0] + 0.5);
			size_t newDatum1 = floor(data[dataNum][1] + 0.5);

			while(newDatum0 >= countTable.size()) {
				if(countTable.size() > 0)
					countTable.push_back(vector<double> (countTable[0].size()));
				else
					countTable.push_back(vector<double>());
			}
			if(newDatum1 >= countTable[0].size())
				for(size_t rowNum = 0; rowNum < countTable.size(); rowNum++)
					countTable[rowNum].resize(newDatum1 + 1);

			countTable[newDatum0][newDatum1]++;
		}

		for(size_t rowNum = 0; rowNum < countTable.size(); rowNum++)
			for(size_t colNum = 0; colNum < countTable[0].size(); colNum++)
				ll += countTable[rowNum][colNum] * (log(countTable[rowNum][colNum] + 1.0) - log((double) data.size() + (double) countTable.size() * (double) countTable[0].size()));
	}

	return(ll);
}
