#ifndef AIS3DTOOLS_PARAMETERS_MANAGER_H
#define AIS3DTOOLS_PARAMETERS_MANAGER_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

namespace Ais3dTools {

class ParametersManager {
  public:
    ////////// CONSTRUCTOR //////////
    ParametersManager(const std::string& filename);
    
    ////////// PUBLIC MEMBER FUNCTIONS //////////
    //! Clears all the saved parameters
    void clear();
    
    //! Gets the current value of a parameter. Returns true if the value was changed
    template <typename ParameterType>
    bool getParameterValue(const std::string& category, const std::string& parameterName, ParameterType& parameter) const;
    
    /** Add/changes the current parameter. Returns true if it either changed or added the parameter
      * \param onlyIfNew controls if the value should also be changed if an entry already exists */
    template <typename ParameterType>
    bool addParameterValue(const std::string& category, const std::string& parameterName, const ParameterType parameter, bool onlyIfNew=false);
    
    //! Checks if the parameters in the file changed
    void updateParameters();
    
    //! Writes the current parameter set to a file
    void writeToFile(const std::string& filename) const;
    
  protected:
    ////////// PROTECTED MEMBER FUNCTIONS //////////
    void readParametersFile();
    
    
    ////////// PROTECTED MEMBER VARIABLES //////////
    typedef std::map<std::string, std::string> ParameterEntries;
    typedef std::map<std::string, ParameterEntries> ParameterEntriesPerCategory;
    ParameterEntriesPerCategory parameters_;
    std::string filename_;
    time_t last_file_modification_time_;
};

} // namespace end

#include "parameters_manager.hpp"

#endif // AIS3DTOOLS_PARAMETERS_MANAGER_H
