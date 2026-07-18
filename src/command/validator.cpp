#include "command/validator.hpp"
#include "core/constants.hpp"

namespace command {

void Validator::validateFileName(
    const std::string& filename)
{
    if (filename.empty()) {
        throw std::runtime_error("Empty filename");
    }
}

void Validator::validateFieldSize(
    int width,
    int height)
{
    validateRange(width,
                        1,
                        std::numeric_limits<int>::max(),
                        "fieldWidth");

    validateRange(height,
                        1,
                        std::numeric_limits<int>::max(),
                        "fieldHeight");
}

void Validator::validateGaussian(
    const DispatcherParams& params)
{
    validateRange(params.centerX,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "centerX");

    validateRange(params.centerY,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "centerY");

    validateRange(params.sigmaX,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sigmaX");

    validateRange(params.sigmaY,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sigmaY");

    validateRange(params.height,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "height");

   
}

void Validator::validateAutoGaussian(
    const DispatcherParams& params)
{
    validateRange(params.xmin,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "xmin");

    validateRange(params.xmax,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "xmax");
                        
   if (params.xmin > params.xmax)
       {
          throw std::runtime_error("xmin must not exceed xmax");
           
       }     
                      
    validateRange(params.ymin,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "ymin");                     

    validateRange(params.ymax,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "ymax");
    
   if (params.ymin > params.ymax)
       {
           throw std::runtime_error("ymin must not exceed ymax");
           
       } 
                           
    validateRange(params.sx_min,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sx_min");
                        
    validateRange(params.sx_max,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sx_max");
         
   if (params.sx_min > params.sx_max)
       {
           throw std::runtime_error("sx_min must not exceed sx_max");
           
       } 
                    
    validateRange(params.sy_min,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sy_min");

    validateRange(params.sy_max,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sy_max");
    
   if (params.sy_min > params.sy_max)
       {
           throw std::runtime_error("sy_min must not exceed sy_max");
           
       } 
                         
    validateRange(params.h_min,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "h_min");

    validateRange(params.h_max,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "h_max");
     
   if (params.h_min > params.h_max)
       {
           throw std::runtime_error("h_min must not exceed h_max");
           
       } 
                        
    validateRange(params.count_min,
                        static_cast<std::size_t>(1),
                        std::numeric_limits<std::size_t>::max(),
                        "count_min");
    
    validateRange(params.count_max,
                        static_cast<std::size_t>(1),
                        std::numeric_limits<std::size_t>::max(),
                        "count_max");                    
    
   if (params.count_min > params.count_max)
       {
           throw std::runtime_error("count_min must not exceed count_max");
           
       }
    
   
}

void Validator::validateGaussGrid(
    int width,
    int height,
    int cellSize)
{
    validateRange(cellSize, 1,
                        std::min(width, height),
                        "gaussCellSize");
                        
    if (width % cellSize != 0 || height % cellSize != 0) {
        throw std::runtime_error("Gauss cellSize must evenly divide both field width and field height.");
    }
}

void Validator::validateBinaryParameters(
    int threshold)
{
    validateRange(threshold,
                         0,
                         static_cast<int>(core::MID_GRAY),
                         "threshold");
}

void Validator::validateBmpWriteMode(
    const std::string& filename,
    const std::string& mode)
{
    validateFileName(filename);
        
    
    if (mode != "Full" &&
        mode != "Binary")
    {    

        throw std::runtime_error("Unknown BMP mode: '" +
            mode +
            "'. Allowed values: Full, Binary");
    }
}

void Validator::validateWaveNoiseSize(
    std::size_t noisy)
{
     validateRange(noisy, static_cast<std::size_t>(0),
                         std::numeric_limits<std::size_t>::max(),
                         "wave noisy");
}

void Validator::validateKMeans(
    std::size_t k)
{
     validateRange(k, static_cast<std::size_t>(1),
                         std::numeric_limits<std::size_t>::max(),
                         "k");
}

void Validator::validateKernelSize(
    std::size_t k,
    std::size_t kernelSize)
{
    validateRange(kernelSize,
                        static_cast<std::size_t>(1),
                        std::numeric_limits<std::size_t>::max(),
                        "kernelSize");

    validateKMeans(k);

   
}

void Validator::validateGrid(
    int width,
    int height,
    int cellSize)
{
    validateRange(cellSize, 1,
                        std::min(width, height),
                        "cellSize");
                        
    if (width % cellSize != 0 || height % cellSize != 0) {
        throw std::runtime_error("Grid cellSize must evenly divide both field width and field height.");
    }
}

void Validator::validateGridNoiseSize(
    std::size_t noisy)
{
     validateRange(noisy, static_cast<std::size_t>(0),
                         std::numeric_limits<std::size_t>::max(),
                         "grid noisy");
}

void Validator::validateNavigationParameters(
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle)
{
    validateRange(vehicleRadius,
                        1,
                        std::numeric_limits<int>::max(),
                        "vehicleRadius");

    validateRange(maxSideAngle,
                        0.0,
                        90.0,
                        "maxSideAngle");

    validateRange(maxUpDownAngle,
                        0.0,
                        90.0,
                        "maxUpDownAngle");

   
}

void Validator::validateConnectParameters(
    const DispatcherParams& params,
    const std::string& mode)
{
    validateConnectParameters(params);
                        
    if (mode != "Nearest" &&
        mode != "NearestK" &&
        mode != "All")
    {
        throw std::runtime_error("Connect mode must be Nearest, NearestK or All");        
    }

    if (mode == "NearestK")
    {
        validateRange(
                params.nearestVerticesCount,
                static_cast<std::size_t>(1),
                std::numeric_limits<std::size_t>::max(),
                "nearestVerticesCount");
    }  
}

void Validator::validateConnectParameters(
    const DispatcherParams& params)
{
    validateRange(params.startPixelX,
                        0,
                        params.fieldWidth - 1,
                        "startPixelX");

    validateRange(params.startPixelY,
                        0,
                        params.fieldHeight - 1,
                        "startPixelY");

    validateRange(params.goalPixelX,
                        0,
                        params.fieldWidth - 1,
                        "goalPixelX");

    validateRange(params.goalPixelY,
                        0,
                        params.fieldHeight - 1,
                        "goalPixelY");
   
}

void Validator::validateRRT(
    const DispatcherParams& params)
{
    validateRange(params.rebuildSize,
                        static_cast<std::size_t>(0),
                        std::numeric_limits<std::size_t>::max(),
                        "rebuildSize");
                        
    validateRange(params.maxIterations,
                        static_cast<std::size_t>(1),
                        std::numeric_limits<std::size_t>::max(),
                        "maxIterations");

    validateRange(params.startWorldX,
                        0.0,
                        static_cast<double>(params.fieldWidth) - 1,
                        "startWorldX");   
                      
    validateRange(params.startWorldY,
                        0.0,
                        static_cast<double>(params.fieldHeight) - 1,
                        "startWorldY");                     

    validateRange(params.goalWorldX,
                        0.0,
                        static_cast<double>(params.fieldWidth) - 1,
                        "goalWorldX");
                           
    validateRange(params.goalWorldY,
                        0.0,
                        static_cast<double>(params.fieldHeight) - 1,
                        "goalWorldY");
                        
    validateRange(params.vehicleRadiusWorld,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::min(params.fieldWidth / 2.0, params.fieldHeight / 2.0)),
                        "vehicleRadiusWorld");
                    
    validateRange(params.heightThresholdWorld,
                        0.0,
                        static_cast<double>(core::MID_GRAY),
                        "heightThresholdWorld");

    validateRange(params.maxSideAngle,
                        0.0,
                        90.0,
                        "maxSideAngle");

    validateRange(params.maxUpDownAngle,
                        0.0,
                        90.0,
                        "maxUpDownAngle");

    validateRange(params.interpEdge,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpEdge");
                        
    validateRange(params.interpCollision,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpCollision");

    validateRange(params.interpAngle,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpAngle");
    
    validateRange(params.step,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "step");

    validateRange(params.goalRadius,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "goalRadius");
    
    validateRange(params.goalBias,
                        0.0,
                        1.0,
                        "goalBias");

   
}

void Validator::validateRRTStar(
    const DispatcherParams& params)
{
    validateRRT(params);
    
    validateRange(params.maxFindRadius,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "maxFindRadius");

    validateRange(params.gammaConstant,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "gammaConstant");
}

void Validator::validateSplineDiscrete(
    const DispatcherParams& params)
{
    validateSplineContinuous(params.samplesPerSegment);
    
    validateRange(params.vehicleRadiusWorld,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::min(params.fieldWidth / 2.0, params.fieldHeight / 2.0)),
                        "vehicleRadiusWorld");
                    
    validateRange(params.heightThresholdWorld,
                        0.0,
                        static_cast<double>(core::MID_GRAY),
                        "heightThresholdWorld");

    validateRange(params.maxSideAngle,
                        0.0,
                        90.0,
                        "maxSideAngle");

    validateRange(params.maxUpDownAngle,
                        0.0,
                        90.0,
                        "maxUpDownAngle");

    validateRange(params.interpEdge,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpEdge");
                        
    validateRange(params.interpCollision,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpCollision");

    validateRange(params.interpAngle,
                        static_cast<double>(core::EPSILON),
                        static_cast<double>(std::max(params.fieldWidth, params.fieldHeight)),
                        "interpAngle");
}

void Validator::validateSplineContinuous(
        std::size_t samplesPerSegment)
{
    validateRange(samplesPerSegment,
                        static_cast<std::size_t>(1),
                        std::numeric_limits<std::size_t>::max(),
                        "samplesPerSegment");
}

} // namespace command
