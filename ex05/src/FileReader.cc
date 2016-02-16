
#include "FileReader.hh"

template<typename T>
void FileReader::registerIntParameter(const std::string &key, int init)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << init;
	mapper[key] = ss.str();}
   //TODO
}

template<typename T>
void FileReader::registerRealParameter(const std::string &key, real init)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << init;
	mapper[key] = ss.str();	}
   //TODO
}

template<typename T>
void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << init;
	mapper[key] = ss.str(); }
   //TODO
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << in;
	mapper[key] = ss.str();}
   //TODO
}

void FileReader::setParameter(const std::string &key, real in)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << in;
	mapper[key] = ss.str();}
   //TODO
}

void FileReader::setParameter(const std::string &key, int in)
{
	if (mapper.find(key) != mapper.end()) {
	WARN("Value already updated")
	std::stringstream ss;
	ss << in;
	mapper[key] = ss.str();}
   //TODO
}


bool FileReader::readFile(const std::string &name)
{
   //TODO
	
	std::ifstream in(name.c_str());
	std::string key,value,buffer;
	if(in){
		while (getline(in, buffer)) {
		std::stringstream ss(buffer);
		ss >> key ; ss >> value ;
		(this->mapper).insert(std::pair<std::string, std::string>(key,value));
	} // End While getline
	in.close(); 
	return true ; }  // end if 
	else {
	CHECK_MSG(true,"**********File Could not be read or Check Again *********************" );
   return false; } // end else	
} // End function



void FileReader::printParameters() const
{
	//Paramap::const_iterator it;
	for (Paramap::const_iterator it = mapper.begin();it != mapper.end();++it){
		std :: cout << it->first << " " << it->second << std::endl;
	}
   //TODO
}


