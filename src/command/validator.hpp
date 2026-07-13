#pragma once

#include <string>
#include <limits>
#include "core/logger.hpp"
#include "command/dispatcher_params.hpp"

namespace command {

class Validator {
private:
    template<typename T>
    static void validateRange(const T& value,
                              const T& minValue,
                              const T& maxValue,
                              const std::string& parameterName)
    {
        if (value < minValue || value > maxValue)
        {
            throw std::runtime_error(
                "Parameter '" + parameterName +
                "' must be in range [" +
                std::to_string(minValue) + ", " +
                std::to_string(maxValue) +
                "], got " +
                std::to_string(value)
            );
        }
    }
    
public:
    static void validateFileName(
        const std::string& filename);
        
    static void validateFieldSize(
        int width,
        int height);
        
    static void validateGaussian(
        const DispatcherParams& params);
        
    static void validateAutoGaussian(
        const DispatcherParams& params);
    
    static void validateBinaryParameters(
        int threshold);
    
    static void validateBmpWriteMode(
        const std::string& filename,
        const std::string& mode);
        
    static void validateWaveNoiseSize(
        std::size_t noisy);

    
    static void validateKMeans(
        std::size_t k);
    
    static void validateKernelSize(
        std::size_t k,
        std::size_t kernelSize);
        
    static void validateGrid(
        int width,
        int height,
        int gridWidth);
        
    static void validateGridNoiseSize(
        std::size_t noisy);
        
    static void validateNavigationParameters(
        int vehicleRadius,
        double maxSideAngle,
        double maxUpDownAngle);   

    static void validateConnectParameters(
        const DispatcherParams& params,
        const std::string& mode);
        
    static void validateConnectParameters(
        const DispatcherParams& params);
        
    static void validateRRT(
        const DispatcherParams& params);
    
    static void validateRRTStar(
        const DispatcherParams& params);
    
    static void validateSplineDiscrete(
        const DispatcherParams& params);
        
    static void validateSplineContinuous(
        std::size_t samplesPerSegment);
};

}
