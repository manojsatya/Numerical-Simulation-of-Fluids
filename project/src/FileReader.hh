#ifndef FILEREADER_HH
#define FILEREADER_HH

#include <string>
#include <map>
#include "Types.hh"
#include "Value.hh"
#include <algorithm>

//*******************************************************************************************************************
/*! Class for reading configuration from a file
*
* Configuration File Syntax:
*   - everything between a '#' character and the beginning of the next line is ignored
*   - lines can be empty
*   - line contain parameter key followed by white-spaces followed by their value
*
*  All possible keys (and the datatype of the value) have to be registered first:
*   - for example usage have a look at the FileReaderTest
*
*
*
*  This Skeleton is just a suggestion. If you are familiar with template programming
*  a possible alternative would be a version that does not need registering of all keys
*  and has a getter function like the following:
*      template<typename ExpectedType>
*      ExpectedType getValue( const std::string & key );
*/
//*******************************************************************************************************************
class FileReader
{
public:

	//register a new parameter with name key and initial int value
	void registerIntParameter( const std::string & key, int init = 0 );

	//register a new parameter with name key and initial double value
	void registerRealParameter( const std::string & key, real init = 0 );

	//register a new parameter with name key and initial string value
	void registerStringParameter( const std::string & key, const std::string & init = "" );

	//set a value for the key string with value in
	void setParameter( const std::string & key, const std::string & in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, real in );

	//set a value for the key string with value in
	void setParameter( const std::string & key, int in );

	// get the int value of key 
	int getIntParameter( const std::string & key )  ;

	// get the double value of key 
	real getRealParameter( const std::string & key ) ;

	// get the string value of key 
	 std::string getStringParameter( const std::string & key ) ;

	Value getParameter(const std::string &key);

	//try to read all registered parameters from file name
	bool readFile( const std::string & name );

	bool checkParameter(const std::string &key);

	//print out all parameters to std:out
	void printParameters() const;

private:
	std::map<std::string,int> IntParameter;
	std::map<std::string,real> RealParameter;
	std::map<std::string,std::string> StringParameter;

};







#endif //FILEREADER_HH

