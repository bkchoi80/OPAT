#ifndef MYUTIL_HPP_
#define MYUTIL_HPP_

#include <vector>
#include <iostream>


namespace myutil{
	bool readDelim(const char* const, std::vector< std::vector<double> >&, bool vervose = false);
}


#endif /*MYUTIL_HPP_*/
