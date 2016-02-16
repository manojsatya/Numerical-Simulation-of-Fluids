#ifndef FILEREADER_HH
#define FILEREADER_HH

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>


#include "Debug.hh"
#include "Types.hh"


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

	// get any value of any type of key using template function 

	template<typename T>
	inline T getParameter(const std::string &key) const ;
/*
	// get the int value of key 
	inline int getIntParameter( const std::string & key ) const;

	// get the double value of key 
	inline real getRealParameter( const std::string & key ) const;

	// get the string value of key 
	inline std::string getStringParameter( const std::string & key ) const; */

	//try to read all registered parameters from file name
	bool readFile( const std::string & name );

	//print out all parameters to std:out
	void printParameters() const;
	
	// get the string input boundary parameters from input file 	
	inline std::string getStringParameterBoundary(const std::string &key) const;

	//get the real boundary paramater from input file 
	inline real getrealParameterBoundary(const std::string &key) const;


private:
	typedef	std::map <std::string,std::string> Paramap;
	Paramap mapper; // Object to map objects
};


template<typename T>
inline T FileReader::getParameter(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		
		std::stringstream ss;T value;
		ss << mapper.find(key)->second;
		ss >> value ;
		return value ;
}
	else  {
		
		CHECK_MSG(false, "************Key not found in the file************");
		
		return -1;	}
   //TODO
   
}

inline std::string FileReader::getStringParameterBoundary(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		std::stringstream ss;std::string value;
		ss << (this->mapper).find(key)->second;
		ss >> value ;
		return value ;
}
	else {
		PROGRESS("*****Key not found , Setting default boundary to No-Slip************");}
   		return "NOSLIP" ;
}


inline real FileReader::getrealParameterBoundary(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		std::stringstream ss; real value ;
		ss << (this->mapper).find(key)->second;
		ss >> value ;
		return value ;
}
	else {
		PROGRESS("*****Key not found , Setting default boundary to 0.0 velocity************");}
		return 0.0 ;
}



/*inline int FileReader::getIntParameter(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		std::stringstream ss;int value;
		ss << mapper.find(key)->second;
		ss >> value ;
		return value ;
}
	else {
		CHECK_MSG(false, "************Key not found in the file************"); }
   //TODO
   
}

inline real FileReader::getRealParameter(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		std::stringstream ss; real value ;
		ss << mapper.find(key)->second;
		ss >> value ;
		return value ;
}
	else {
		CHECK_MSG(false, "************Key not found in the file************");}
   //TODO
   //return 0.0;
}

inline std::string FileReader::getStringParameter(const std::string &key) const
{
	if (mapper.find(key) != mapper.end()) {
		std::stringstream ss; std::string value ;
		ss << mapper.find(key)->second;
		ss >> value ;
		return value ;
}
	else {
		CHECK_MSG(false, "************Key not found in the file************");}
   //TODO
//   return "";
}

*/



#endif //FILEREADER_HH

