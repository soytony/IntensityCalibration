#include "parameters_manager.h"

#include <sys/stat.h>
//#include <unistd.h>
//#include <time.h>

////////////////////////////////////////////////////////////

Ais3dTools::ParametersManager::ParametersManager(const std::string& filename) : filename_(filename), last_file_modification_time_(0)
{
}

////////////////////////////////////////////////////////////

void Ais3dTools::ParametersManager::updateParameters() {
  struct stat fileStats;
  stat(filename_.c_str(), &fileStats);
  if (fileStats.st_mtime != last_file_modification_time_) {
    //std::cout << "File changed.\n";
    readParametersFile();
    last_file_modification_time_ = fileStats.st_mtime;
  }
  else {
    //std::cout << "File remains unchanged.\n";
  }
}

////////////////////////////////////////////////////////////

void Ais3dTools::ParametersManager::readParametersFile() {
  std::ifstream file(filename_.c_str());
  if (!file) {
    std::cerr << "Could not read parameters file \""<<filename_<<"\".\n";
    return;
  }
  std::string currentCategory = "global";
  
  while (file.peek() != EOF) {
    std::string tmpLine;
    std::getline(file, tmpLine);
    std::stringstream line(tmpLine);
    std::string token;
    line >> token;
    //std::cout << PVARN(token);
    if (token.empty())
      continue;
    if (token[0]=='[') {
      if (token[token.size()-1] != ']' || token.size()<3) {
        std::cerr << "Error reading category tag \""<<token<<"\".\n";
        continue;
      }
      currentCategory = token.substr(1,token.size()-2);
      //std::cout << PVARN(currentCategory);
      continue;
    }
    
    if (token[0]=='#')  // Line is a comment?
      continue;
    
    std::string parameterValue;
    line >> parameterValue;
    if (parameterValue.empty()) {
      std::cerr << "Error reading value of parameter \""<<token<<"\" in category \""<<currentCategory<<"\".\n";
      continue;
    }
    //std::cout << "Setting value of parameter \""<<token<<"\" in category \""<<currentCategory<<"\" to \""<<parameterValue<<"\".\n";
    parameters_[currentCategory][token] = parameterValue;
  }
  file.close();
  
  //std::cout << "\nCurrent parameters\n------------------\n";
  //for (ParameterEntriesPerCategory::const_iterator it=parameters_.begin(); it!=parameters_.end(); ++it) {
    //std::cout << "["<<it->first<<"]\n";
    //for (ParameterEntries::const_iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
      //std::cout << "  "<<it2->first<<" = "<<it2->second<<"\n";
    //}
  //}
  //std::cout << "\n\n";
}

////////////////////////////////////////////////////////////

void Ais3dTools::ParametersManager::clear() {
  parameters_.clear();
}

////////////////////////////////////////////////////////////

void Ais3dTools::ParametersManager::writeToFile(const std::string& filename) const {
  std::ofstream file(filename.c_str(), std::ios_base::out|std::ios_base::trunc);
  if (!file) {
    std::cerr << "Could not open file \""<<filename<<"\" for writing.\n";
  }
  
  for (ParameterEntriesPerCategory::const_iterator it=parameters_.begin(); it!=parameters_.end(); ++it) {
    file << "["<<it->first<<"]\n";
    for (ParameterEntries::const_iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
      file << "  "<<it2->first<<" "<<it2->second<<"\n";
    }
  }
}
