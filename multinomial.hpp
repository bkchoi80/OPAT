/*
 * multinomial.hpp
 *
 *  Created on: May 25, 2010
 *      Author: bkchoi
 */

#ifndef MULTINOMIAL_HPP_
#define MULTINOMIAL_HPP_

#include <vector>

namespace multinomial
{
	double computeLogLikelihood(const std::vector< std::vector<double> >&);
}

#endif /* MULTINOMIAL_HPP_ */
