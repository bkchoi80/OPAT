#include <limits>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <assert.h>

#include "opat.hpp"
#include "myutil.hpp"


void OPAT::initialize(const std::vector< std::vector<double> > &data)
{
	// Discretize the data with 15% padding at both of the ends.
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
		double minPos = data[0][dimNum];
		double maxPos = minPos;
		for(size_t dataNum = 1; dataNum < data.size(); dataNum++) {
			minPos = (minPos > data[dataNum][dimNum]) ? data[dataNum][dimNum] : minPos;
			maxPos = (maxPos < data[dataNum][dimNum]) ? data[dataNum][dimNum] : maxPos;
		}

		mOffset[dimNum] = minPos - (maxPos - minPos) * 0.15;
		mScale[dimNum] = maxPos - minPos + 2.0 * (maxPos - minPos) * 0.15;
	}

	for(size_t dataNum = 0; dataNum < mData.size(); dataNum++) {
		mData[dataNum].resize(mDimension);
		for(size_t dimNum = 0; dimNum < mDimension; dimNum++)
			mData[dataNum][dimNum] = round(pow(2, MAX_SPLIT_DEPTH) * (data[dataNum][dimNum] - mOffset[dimNum]) / mScale[dimNum]);
	}

	mElRegions.clear();
}

void OPAT::initialize(const std::vector< std::vector<double> > &data, const std::vector<size_t> &varIndex)
{
	for(size_t index = 0; index < mDimension; index++) {
		size_t dimNum = varIndex[index];
		double minPos = data[0][dimNum];
		double maxPos = minPos;
		for(size_t dataNum = 1; dataNum < data.size(); dataNum++) {
			minPos = (minPos > data[dataNum][dimNum]) ? data[dataNum][dimNum] : minPos;
			maxPos = (maxPos < data[dataNum][dimNum]) ? data[dataNum][dimNum] : maxPos;
		}

		mOffset[index] = minPos - (maxPos - minPos) / M_PI;
		mScale[index] = maxPos - minPos + 2.0 * (maxPos - minPos) / M_PI;
	}

	for(size_t dataNum = 0; dataNum < mData.size(); dataNum++) {
		mData[dataNum].resize(mDimension);
		for(size_t index = 0; index < mDimension; index++) {
			size_t dimNum = varIndex[index];
			mData[dataNum][index] = round(pow(2, MAX_SPLIT_DEPTH) * (data[dataNum][dimNum] - mOffset[index]) / mScale[index]);
		}
	}

	mElRegions.clear();
}

void OPAT::build(size_t numLookAheads, double rho)
{
	OPATNode root(mDimension);
	root.mPoints.resize(mData.size());
	for(size_t dataNum = 0; dataNum < mData.size(); dataNum++)
		root.mPoints[dataNum] = dataNum;

	buildSub(root, numLookAheads, rho);
}

void OPAT::buildSub(OPATNode &node, size_t numLookAheads, double rho)
{
	// Fully grow the current tree by numLookAheads levels.
	buildLookAhead(node, numLookAheads, rho);

	// If the partition as one point or the stopping posterior is greater than 0.5 than stop and store it to the ElRegions vector.
	// Otherwise, choose the best split by the posterior lambdas after lookahead.
	if(node.mPoints.size() < 2 || node.mRho > 0.5) {
		node.mLeftChildren.clear();
		node.mRightChildren.clear();
		mElRegions.push_back(node);
	}
	else {
		size_t bestSplitIndex = std::max_element(node.mLambdas.begin(), node.mLambdas.end()) - node.mLambdas.begin();
		node.mLeftChildren.erase(node.mLeftChildren.begin(), node.mLeftChildren.begin() + bestSplitIndex);
		node.mLeftChildren.erase(node.mLeftChildren.begin() + 1, node.mLeftChildren.end());
		node.mRightChildren.erase(node.mRightChildren.begin(), node.mRightChildren.begin() + bestSplitIndex);
		node.mRightChildren.erase(node.mRightChildren.begin() + 1, node.mRightChildren.end());

		node.mPoints.clear();

		node.mLeftChildren[0].mProbMass = node.mProbMass * (node.mLeftChildren[0].mPoints.size() + 0.5) / (node.mLeftChildren[0].mPoints.size() + node.mRightChildren[0].mPoints.size() + 1.0);
		node.mRightChildren[0].mProbMass = node.mProbMass * (node.mRightChildren[0].mPoints.size() + 0.5) / (node.mLeftChildren[0].mPoints.size() + node.mRightChildren[0].mPoints.size() + 1.0);
		buildSub(node.mLeftChildren[0], numLookAheads, rho);
		buildSub(node.mRightChildren[0], numLookAheads, rho);
	}
}

void OPAT::buildLookAhead(OPATNode &node, size_t numLookAheads, double rho)
{
	assert(node.mLeftChildren.size() == node.mRightChildren.size() && node.mRightChildren.size() == node.mLambdas.size());

	if(numLookAheads == 0 || node.mPoints.size() < 2)
		return;

	// If the split has not been computed than assign data into children.
	if(node.mLambdas.size() == 0) {
		for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
			if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
				OPATNode leftChild(mDimension);
				leftChild.mDepths = node.mDepths;
				leftChild.mPosition = node.mPosition;
				leftChild.mDepths[dimNum]++;
				leftChild.mPosition[dimNum] = node.mPosition[dimNum] << 1;

				OPATNode rightChild(mDimension);
				rightChild.mDepths = node.mDepths;
				rightChild.mPosition = node.mPosition;
				rightChild.mDepths[dimNum]++;
				rightChild.mPosition[dimNum] = (node.mPosition[dimNum] << 1) + 1;

				for(size_t dataNum = 0; dataNum < node.mPoints.size(); dataNum++)
					if((mData[node.mPoints[dataNum]][dimNum] >> (MAX_SPLIT_DEPTH - node.mDepths[dimNum] - 1)) & 1u)
						rightChild.mPoints.push_back(node.mPoints[dataNum]);
					else
						leftChild.mPoints.push_back(node.mPoints[dataNum]);

				node.mLeftChildren.push_back(leftChild);
				node.mRightChildren.push_back(rightChild);
			}
		}

		node.mLambdas.resize(mDimension);
	}

	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
		if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
			buildLookAhead(node.mLeftChildren[dimNum], numLookAheads - 1, rho);
			buildLookAhead(node.mRightChildren[dimNum], numLookAheads - 1, rho);
		}
	}

	// Do calculations in the log space for numerical stability.
	double phi = 0;
	std::vector<double> secondTerms(mDimension);
	double maxSecondTerm;
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
		if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
			secondTerms[dimNum] = computeCoeff(node.mLeftChildren[dimNum].mPoints.size(), node.mRightChildren[dimNum].mPoints.size()) + node.mLeftChildren[dimNum].mLogPhi + node.mRightChildren[dimNum].mLogPhi;
			maxSecondTerm = secondTerms[dimNum];
		}
	}
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++)
			if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH)
				maxSecondTerm = (maxSecondTerm < secondTerms[dimNum]) ? secondTerms[dimNum] : maxSecondTerm;
	double sum = 0;
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
		if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
			secondTerms[dimNum] -= maxSecondTerm;
			sum += exp(secondTerms[dimNum]);
		}
	}
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
		if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
			node.mLambdas[dimNum] = exp(secondTerms[dimNum]) / sum;
			phi += exp(secondTerms[dimNum]);
		}
		else
			node.mLambdas[dimNum] = 0.0;
	}
	bool fullSplit = true;
	for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
			if(node.mDepths[dimNum] < MAX_SPLIT_DEPTH) {
				fullSplit = false;
				break;
			}
	}

	if(fullSplit) {
		node.mLogPhi = 0.0;
	}
	else {
		if(maxSecondTerm > 0.0) {
			phi = rho / exp(maxSecondTerm) + phi * (1-rho) / mDimension;
			node.mLogPhi = maxSecondTerm + log(phi);
		}
		else {
			node.mLogPhi = log(rho + exp(maxSecondTerm) * phi * (1-rho) / mDimension);
		}
		node.mRho = rho / exp(node.mLogPhi);
	}
}

inline double OPAT::computeCoeff(size_t leftNum, size_t rightNum)
{
	double logVal = lgamma(0.5 + leftNum) + lgamma(0.5 + rightNum) - lgamma(1 + leftNum + rightNum) - lgamma(M_PI);
	logVal += (leftNum + rightNum) * log(2);
	return(logVal);
}

void OPAT::writePostDensity(std::string filename)
{
	std::ofstream ofs;
	ofs.open(filename.c_str(), std::ios_base::app);

	for(size_t regionNum = 0; regionNum < mElRegions.size(); regionNum++) {
		OPATNode &region = mElRegions[regionNum];

		double mu = 1.0;
		for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
			double lowerBound, upperBound;
			if(region.mDepths[dimNum] == 0) {
				lowerBound = mOffset[dimNum];
				upperBound = lowerBound + mScale[dimNum];
			}
			else {
				lowerBound = mOffset[dimNum] + mScale[dimNum] *  (region.mPosition[dimNum] << (MAX_SPLIT_DEPTH - region.mDepths[dimNum])) / pow(2, MAX_SPLIT_DEPTH);
				upperBound = lowerBound + mScale[dimNum] / pow(2, region.mDepths[dimNum]);
			}
			mu *= upperBound - lowerBound;
			ofs << lowerBound << "\t" << upperBound << "\t";
		}

		ofs << region.mProbMass / mu;
		ofs << std::endl;
	}

	ofs.close();
}

double OPAT::computeConditionalML(size_t var)
{
	assert(var < mDimension);

	double ml = 0.0;
	for(size_t regionNum = 0; regionNum < mElRegions.size(); regionNum++) {
		double depth = 0.0;
		for(size_t dimNum = 0; dimNum < mDimension; dimNum++)
			depth += mElRegions[regionNum].mDepths[dimNum];
		ml += mElRegions[regionNum].mPoints.size() * (log(mElRegions[regionNum].mProbMass) + depth * log(2));
	}

	double marginalML = 0.0;
	for(size_t dataNum = 0; dataNum < mData.size(); dataNum++) {
		std::vector<double> marginals;
		for(size_t regionNum = 0; regionNum < mElRegions.size(); regionNum++) {
			OPATNode &region = mElRegions[regionNum];

			bool covered = true;
			for(size_t dimNum = 0; dimNum < mDimension; dimNum++) {
				if((dimNum != var) && (region.mDepths[dimNum] != 0) && ((mData[dataNum][dimNum] >> (MAX_SPLIT_DEPTH - region.mDepths[dimNum])) != region.mPosition[dimNum])) {
					covered = false;
					break;
				}
			}

			if(covered) {
				marginals.push_back(log(region.mProbMass));
				for(size_t dimNum = 0; dimNum < mDimension; dimNum++)
					if(dimNum != var)
						marginals.back() += log(2) * region.mDepths[dimNum];
			}
		}

		double maxMarginal = *max_element(marginals.begin(), marginals.end());
		double marginalSum = 0.0;
		for(size_t marginalNum = 0; marginalNum < marginals.size(); marginalNum++)
			marginalSum += exp(marginals[marginalNum] - maxMarginal);
		marginalSum = maxMarginal + log(marginalSum);
		marginalML += marginalSum;
	}

	return(ml - marginalML - mData.size() * log(mScale[var]));
}

void OPAT::debug()
{

}
