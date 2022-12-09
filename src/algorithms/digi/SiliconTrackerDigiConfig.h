#ifndef EICRECON_SILICONTRACKERDIGICONFIG_H
#define EICRECON_SILICONTRACKERDIGICONFIG_H

#include <DD4hep/DD4hepUnits.h>

namespace eicrecon {

    struct SiliconTrackerDigiConfig {
        double threshold  = 0;
        double timeResolution = 8*dd4hep::ns;
    };

} // eicrecon

#endif //EICRECON_SILICONTRACKERDIGICONFIG_H
