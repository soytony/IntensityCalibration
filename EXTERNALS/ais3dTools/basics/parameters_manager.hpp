////////////////////////////////////////////////////////////

template <typename ParameterType>
bool Ais3dTools::ParametersManager::getParameterValue(const std::string& category, const std::string& parameterName, ParameterType& parameter) const {
  const ParameterEntriesPerCategory::const_iterator it = parameters_.find(category);
  if (it==parameters_.end()) {
    std::cerr << "Unknown category \""<<category<<"\".\n";
    return false;
  }
  const ParameterEntries::const_iterator it2 = it->second.find(parameterName);
  if (it2==it->second.end()) {
    std::cerr << "Unknown parameter \""<<parameterName<<"\" in category \""<<category<<"\".\n";
    return false;
  }
  ParameterType newParamValue;
  std::stringstream ss(it2->second);
  ss >> newParamValue;
  //std::cout << "Parameter \""<<category<<"_"<<parameterName<<"\" was set to "<<parameter<<" from string value \""<<it2->second<<"\".\n";
  if (newParamValue == parameter) {
    //std::cout << "Parameter \""<<category<<"_"<<parameterName<<"\" remains unchanged at value "<<parameter<<".\n";
    return false;
  }
  //std::cout << "Parameter \""<<category<<"_"<<parameterName<<"\" was changed from "<<parameter<<" to "<<newParamValue<<".\n";
  parameter = newParamValue;
  return true;
}

template <typename ParameterType>
bool Ais3dTools::ParametersManager::addParameterValue(const std::string& category, const std::string& parameterName,
                                                      const ParameterType parameter, bool onlyIfNew)
{
  ParameterEntries& entriesInCategory = parameters_[category];
  const ParameterEntries::const_iterator it = entriesInCategory.find(parameterName);
  if (it==entriesInCategory.end()) {
    std::stringstream tmpStringStream;
    tmpStringStream << parameter;
    entriesInCategory[parameterName] = tmpStringStream.str();
    return true;
  }
  
  if (onlyIfNew)
    return false;
  
  ParameterType oldParamValue;
  std::stringstream tmpStringStream(it->second);
  tmpStringStream >> oldParamValue;
  if (oldParamValue == parameter)
    return false;
  
  //std::cout << "Parameter \""<<category<<"_"<<parameterName<<"\" was changed from "<<oldParamValue<<" to "<<parameter<<".\n";
  
  tmpStringStream.str("");
  tmpStringStream << parameter;
  entriesInCategory[parameterName] = tmpStringStream.str();
  return true;
}
