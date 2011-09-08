#include <iostream>
#include <boost/spirit.hpp>

#include "myutil.hpp"


namespace myutil {

	using namespace boost::spirit;
	using namespace std;
	typedef file_iterator<char> f_iterator_t;
	typedef scanner<f_iterator_t> f_scanner_t;
	typedef rule<f_scanner_t> f_rule_t;

	struct TempData
	{
		vector< vector<double> >& data;
		vector<double> datapoint;
		unsigned long n;
		unsigned long counter;
	    TempData(vector< vector<double> >& data_)
	    : data(data_), datapoint(0), n(0), counter(0) {};
	};

	struct update
	{
	    update(TempData& temp_data_)
	    : temp_data(temp_data_) {}

	    void operator()(f_iterator_t first, f_iterator_t last) const
	    {
	   		if( temp_data.counter != temp_data.datapoint.size() ){
	   			cout << "Different number of entries. Parse failed!\n";
	        	exit(-1);
	        }
	   		temp_data.data.push_back(temp_data.datapoint);
	   		temp_data.counter = 0;
	   		temp_data.n++;
	    }

	    TempData& temp_data;
	};

	struct update_first
	{
	    update_first(TempData& temp_data_)
	    : temp_data(temp_data_) {}

	    void operator()(f_iterator_t first, f_iterator_t last) const
	    {
	   		temp_data.data.push_back(temp_data.datapoint);
	   		temp_data.n++;
	    }

	    TempData& temp_data;
	};

	struct assign_entry
	{
	    assign_entry(TempData& temp_data_)
	    : temp_data(temp_data_) {}

	    void operator()(double entry) const
	    {
	        temp_data.datapoint[temp_data.counter] = entry;
	        temp_data.counter++;
	        if(temp_data.counter > temp_data.datapoint.size()){
	        	cout << "Different number of entries. Parse failed!\n";
	        	exit(-1);
	        }
	    }

	  	TempData& temp_data;
	};

	bool readDelim(const char* filename, vector< vector<double> >& data, bool verbose){
	    data.clear();

	    f_iterator_t first(filename);
	    if (!first)
	    {
	        std::cout << "Unable to open file!\n";
	        return -1;
	    }
	    f_iterator_t last = first.make_end();
		TempData tempdata(data);

	    f_rule_t first_line =
	    	list_p.direct(real_p[push_back_a(tempdata.datapoint)], chset<>(" ,\t"))
	    	>> eol_p >> eps_p[update_first(tempdata)];
	    f_rule_t line =
	    	list_p.direct(real_p[assign_entry(tempdata)], chset<>(" ,\t"))
	    	>> eol_p >> eps_p[update(tempdata)];

	    parse_info <f_iterator_t> info = parse(first, last, first_line >> *line >> *space_p);


	    if (info.full){
	    	if(verbose) {
	    		cout << "Parse succeeded!\n";
	    		cout << tempdata.datapoint.size() << " variable data of size " << tempdata.n << "\n";
	    	}
	    }
	    else{
	        std::cout << "Parse failed!\n";
	        return -1;
	    }

	    return 0;
	}

}
