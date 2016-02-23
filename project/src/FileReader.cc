
#include "FileReader.hh"
#include <fstream>
#include "Debug.hh"
#include <sstream>
#include <limits>
#include <iostream>

//TODO implement a method for case insensitive reading of string

void FileReader::registerIntParameter(const std::string &key, int init)
{
   IntParameter[key] = 0;
   setParameter(key,init);
}

void FileReader::registerRealParameter(const std::string &key, real init)
{
   RealParameter[key] = 0.0;
   setParameter(key,init);
}

void FileReader::registerStringParameter(const std::string &key, const std::string &init)
{
   StringParameter[key] = " ";
   setParameter(key,init);
}

void FileReader::setParameter(const std::string &key, const std::string &in)
{
	 StringParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, real in)
{
	 RealParameter[key] = in;
}

void FileReader::setParameter(const std::string &key, int in)
{
	IntParameter[key] = in;
}

Value FileReader::getParameter(const std::string &key)
{
	 std::stringstream ss,msg;
	auto search_string = StringParameter.find(key);
	auto search_int    = IntParameter.find(key);
	auto search_real   =  RealParameter.find(key);
	if(search_string != StringParameter.end())
	{
		StringParameter[key].erase( std::remove(StringParameter[key].begin(), StringParameter[key].end(), ' '), StringParameter[key].end() );
		ss<<StringParameter[key];
		return(ss.str());
	}
	else if(search_int != IntParameter.end() )
	{
		ss<<IntParameter[key];
		return(ss.str());
	}
	else if(search_real != RealParameter.end() )
	{
		ss<<RealParameter[key];
		return(ss.str());;
	}
    else
    {
    	msg<<"Parameter "<<key<<" not found";
    	WARN(msg.str());
	   return (ss.str());
    }
}
//returns true if key present else false
bool FileReader::checkParameter(const std::string &key)
{
	auto  search_int    = IntParameter.find(key);
	auto  search_real   = RealParameter.find(key);
	auto search_string  = StringParameter.find(key);

	 if(search_int == IntParameter.end() && search_real == RealParameter.end() && search_string == StringParameter.end())
	 {
		 return false;
	 }
	 else
		 return true;

}
int FileReader::getIntParameter(const std::string &key)
{
 auto  search_int    = IntParameter.find(key);
 if(search_int != IntParameter.end())
 return IntParameter[key];
 else
 {
	std::stringstream msg;
	msg<<"Parameter " <<key<< " not found, Setting to invalid int 0";
	WARN(msg.str());
 return 0;
 }
}

real FileReader::getRealParameter(const std::string &key)
{
auto  search_int    = RealParameter.find(key);
if(search_int != RealParameter.end())
{
return RealParameter[key];

}
else
{
std::stringstream msg;
msg<<"Parameter " <<key<< " not found, Setting to invalid real 0.0";
WARN(msg.str());
return 0.0;
}

}

std::string FileReader::getStringParameter(const std::string &key)
{
auto  search_int    = StringParameter.find(key);
if(search_int != StringParameter.end())
{
  StringParameter[key].erase( std::remove(StringParameter[key].begin(), StringParameter[key].end(), ' '), StringParameter[key].end() );
  return StringParameter[key];

}
else
{
std::stringstream msg;
msg<<"Parameter " <<key<< " not found, Setting to invalid string  ";
WARN(msg.str());
return "";
}
}



bool FileReader::readFile(const std::string &name)
{
   PRG_LEVEL("Reading Input file");
   std::fstream file (name,std::fstream::in );
   std::stringstream ss;
   ss << "Input File "<<name<<" cannot be opened";
   std::string err_msg = ss.str();
   ASSERT_MSG(file, err_msg);

   bool flag = false;
   int int_inp;
   real real_inp;
   std::string string_inp, line, key;

   std::size_t pos,char_pos,comment_pos,tab_pos,space_pos;

   bool blank_line;

   if(file.good())
	   flag =true;

   while(file.good())
   {
	  getline(file,line);
	  blank_line = false;
	  char_pos = 0;
	  size_t i = 0;
	  while(!isalpha(line[i])&& i<line.size())
	       ++i;
	  char_pos = i;
	  if(char_pos == line.size())
	  {
		  blank_line = true;
	  }

	  comment_pos = line.size();
	  comment_pos = line.find('#');

	  line = line.substr(char_pos);

	  if((comment_pos > char_pos) && !blank_line)
	  {
	   tab_pos = line.find('\t');
	   space_pos = line.find(' ');
	   pos = std::min(space_pos,tab_pos);
	   key = line.substr(0,pos);
	   size_t i = pos+1;
	   while(line[i]==' '||line[i]=='\t')
		  ++i;
	   pos = i;
	   std::stringstream iss ;
	   iss << line.substr(pos,comment_pos-pos);

	   if(iss >> real_inp)
	   	  {
		  if(real_inp == static_cast<int>(real_inp))
		  {
			  int_inp = real_inp;
			  setParameter(key,int_inp);

		  }
		  else
		  {
			  setParameter(key,real_inp);

	   	  }
	   	  }
	   else
		  {
		   setParameter(key,line.substr(pos,comment_pos-pos-1));
		  }
	  }
     }
   return flag;
 }

void FileReader::printParameters() const
{
   for(auto iter = IntParameter.begin(); iter != IntParameter.end(); ++iter)
   {
	   std::cout<< iter ->first<<"\t"<<iter ->second<<std::endl;
   }
   for(auto iter = RealParameter.begin(); iter != RealParameter.end(); ++iter)
   {
	   std::cout<< iter ->first<<"\t"<<iter ->second<<std::endl;
   }
   for(auto iter = StringParameter.begin(); iter != StringParameter.end(); ++iter)
   {
	   std::cout<< iter ->first<<"\t"<<iter ->second<<std::endl;
   }
}


