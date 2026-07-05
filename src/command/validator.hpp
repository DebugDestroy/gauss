#pragma once

#include <string>
#include <limits>
#include "core/logger.hpp"
#include "command/dispatcher_params.hpp"

namespace command {

class Validator {
private:
    template<typename T>
    static bool validateRange(const T& value,
                              const T& minValue,
                              const T& maxValue,
                              const std::string& parameterName,
                              core::Logger& logger)
    {
        if (value < minValue || value > maxValue)
        {
            logger.error(
                "Parameter '" + parameterName +
                "' must be in range [" +
                std::to_string(minValue) + ", " +
                std::to_string(maxValue) + "]");

            return false;
        }

        return true;
    }
    
public:
    static bool validateFileName(
        const std::string& filename,
        core::Logger& logger);
        
    static bool validateFieldSize(
        int width,
        int height,
        core::Logger& logger);
        
    static bool validateGaussian(
        const DispatcherParams& params,
        core::Logger& logger);
        
    static bool validateAutoGaussian(
        const DispatcherParams& params,
        const std::string& mode,
        core::Logger& logger);
    
    static bool validateBinaryParameters(
        int threshold,
        const std::string& mode,
        core::Logger& logger);
    
    static bool validateBmpWriteMode(
        const std::string& filename,
        const std::string& mode,
        core::Logger& logger);
        
    static bool validateWaveNoiseSize(
        int noisy,
        core::Logger& logger);

    
    static bool validateKMeans(
        int k,
        core::Logger& logger);
    
    static bool validateKernelSize(
        int k,
        int kernelSize,
        core::Logger& logger);
        
    static bool validateGrid(
        int width,
        int height,
        int gridWidth,
        core::Logger& logger);
        
    static bool validateGridNoiseSize(
        int noisy,
        core::Logger& logger);
        
    static bool validateNavigationParameters(
        int vehicleRadius,
        double maxSideAngle,
        double maxUpDownAngle,
        core::Logger& logger);   

    static bool validateConnectParameters(
        const DispatcherParams& params,
        const std::string& mode,
        core::Logger& logger);
        
    static bool validateConnectParameters(
        const DispatcherParams& params,
        core::Logger& logger);
};

}
