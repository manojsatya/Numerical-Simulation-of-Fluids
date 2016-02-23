#include <sstream>
#include "Debug.hh"
struct Value
{
std::string value;
Value(std::string value_):value(value_){};

template<typename T> operator T() const
 {
	std::stringstream ss(value);
	std::stringstream msg_ss;
	T out_value ;
	if(ss>>out_value)
		return out_value;
	else
	{
		msg_ss<<"THE RETURN TYPE OF THE KEY VALUE DOES NOT MATCH, RIGHT SIDE OR NO SUCH PARAMETER";
		ERROR(msg_ss.str());
		return 0;
	}
 }

};
