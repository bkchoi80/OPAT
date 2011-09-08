#ifndef OPAT_HPP_
#define OPAT_HPP_

#include <vector>
#include <fstream>
#include <string>

static const size_t MAX_SPLIT_DEPTH = 32;


class OPATNode
{
	friend class OPAT;

public:
	OPATNode(size_t dimension) : mDimension(dimension), mPosition(dimension), mDepths(dimension), mLambdas(0), mPoints(0), mLogPhi(0.0), mRho(1.0), mProbMass(1.0) {};

private:
	std::vector<OPATNode> mLeftChildren;
	std::vector<OPATNode> mRightChildren;
	size_t mDimension;
	std::vector<size_t> mPosition;
	std::vector<size_t> mDepths;
	std::vector<double> mLambdas;
	std::vector<size_t> mPoints;
	double mLogPhi;
	double mRho;
	double mProbMass;
};

class OPAT
{
public:
	OPAT(const std::vector< std::vector<double> > &data) : mDimension(data[0].size()), mOffset(mDimension), mScale(mDimension), mData(data.size()) { initialize(data); };
	void build(size_t, double);
	void writePostDensity(const std::string filename);
	double computeConditionalML(size_t var);
	inline double computeML();
	double computeMI();
	double computeEntropy();
	void debug();

private:
	size_t mDimension;
	std::vector<double> mOffset;
	std::vector<double> mScale;
	std::vector< std::vector<size_t> > mData;
	std::vector<OPATNode> mElRegions;

	void initialize(const std::vector< std::vector<double> > &data);
	inline double computeCoeff(size_t leftNum, size_t rightNum);
	void buildSub(OPATNode &node, size_t numLookAheads, double rho);
	void buildLookAhead(OPATNode &node, size_t numLookAheads, double rho);
};


#endif /* OPAT_HPP_ */
