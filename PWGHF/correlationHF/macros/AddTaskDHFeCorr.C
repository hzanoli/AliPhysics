#include "AliAnalysisTaskDHFeCorr.h"
#include <string>

AliAnalysisTaskDHFeCorr *AddTaskDHFeCorr(std::string name, std::string config_file, Int_t trigger = AliVEvent::kINT7) {
    return new AliAnalysisTaskDHFeCorr(name.c_str());
}
