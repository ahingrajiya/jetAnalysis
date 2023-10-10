#ifndef BASEHISTOMANAGER_H
#define BASEHISTOMANAGER_H

// ROOT headers
#include "TList.h"

//________________
class BaseHistoManager {
  public:
    /// @brief Constructor
    BaseHistoManager() { /* empty */ }
    /// @brief Destructor
    virtual ~BaseHistoManager() { /* empty */}
    
    /// @brief Return list of objects (histograms, profiles, graphs, etc...)
    virtual TList *getOutputList() = 0;

  ClassDef(BaseHistoManager, 0)
};

#endif // #define BASEHISTOMANAGER_H