#include "command/validator.hpp"
#include "core/constants.hpp"

namespace command {

bool Validator::validateFileName(
    const std::string& filename,
    core::Logger& logger)
{
    if (filename.empty()) {
        logger.error("Filename cannot be empty");
        return false;
    }

    return true;
}

bool Validator::validateFieldSize(
    int width,
    int height,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(width,
                        1,
                        std::numeric_limits<int>::max(),
                        "fieldWidth",
                        logger);

    ok &= validateRange(height,
                        1,
                        std::numeric_limits<int>::max(),
                        "fieldHeight",
                        logger);

    return ok;
}

bool Validator::validateGaussian(
    const DispatcherParams& params,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(params.centerX,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "centerX",
                        logger);

    ok &= validateRange(params.centerY,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "centerY",
                        logger);

    ok &= validateRange(params.sigmaX,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sigmaX",
                        logger);

    ok &= validateRange(params.sigmaY,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sigmaY",
                        logger);

    ok &= validateRange(params.height,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "height",
                        logger);

    return ok;
}

bool Validator::validateAutoGaussian(
    const DispatcherParams& params,
    const std::string& mode,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(params.xmin,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "xmin",
                        logger);

    ok &= validateRange(params.xmax,
                        0.0,
                        static_cast<double>(params.fieldWidth),
                        "xmax",
                        logger);
                        
   if (params.xmin > params.xmax)
       {
           logger.error("xmin must not exceed xmax");
           ok = false;
       }     
                      
    ok &= validateRange(params.ymin,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "ymin",
                        logger);                     

    ok &= validateRange(params.ymax,
                        0.0,
                        static_cast<double>(params.fieldHeight),
                        "ymax",
                        logger);
    
   if (params.ymin > params.ymax)
       {
           logger.error("ymin must not exceed ymax");
           ok = false;
       } 
                           
    ok &= validateRange(params.sx_min,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sx_min",
                        logger);
                        
    ok &= validateRange(params.sx_max,
                        static_cast<double>(core::EPSILON),
                        params.fieldWidth / 2.0,
                        "sx_max",
                        logger);
         
   if (params.sx_min > params.sx_max)
       {
           logger.error("sx_min must not exceed sx_max");
           ok = false;
       } 
                    
    ok &= validateRange(params.sy_min,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sy_min",
                        logger);

    ok &= validateRange(params.sy_max,
                        static_cast<double>(core::EPSILON),
                        params.fieldHeight / 2.0,
                        "sy_max",
                        logger);
    
   if (params.sy_min > params.sy_max)
       {
           logger.error("sy_min must not exceed sy_max");
           ok = false;
       } 
                         
    ok &= validateRange(params.h_min,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "h_min",
                        logger);

    ok &= validateRange(params.h_max,
                        static_cast<double>((-1)*core::MID_GRAY),
                        static_cast<double>(core::MID_GRAY),
                        "h_max",
                        logger);
     
   if (params.h_min > params.h_max)
       {
           logger.error("h_min must not exceed h_max");
           ok = false;
       } 
                        
    ok &= validateRange(params.count_min,
                        1,
                        std::numeric_limits<int>::max(),
                        "count_min",
                        logger);
    
    ok &= validateRange(params.count_max,
                        1,
                        std::numeric_limits<int>::max(),
                        "count_max",
                        logger);                    
    
   if (params.count_min > params.count_max)
       {
           logger.error("count_min must not exceed count_max");
           ok = false;
       } 
   
   if (mode != "Random" &&
       mode != "Fixed")
    {
        logger.error("GAuto mode must be Random, Fixed");

        ok = false;
    }

    if (mode == "Fixed")
    {
        if (!validateRange(
                params.seedGAuto,
                0u,
                std::numeric_limits<uint32_t>::max(),
                "seedGAuto",
                logger))
        {
            ok = false;
        }
    }
    
    return ok;
}

bool Validator::validateBinaryParameters(
    int threshold,
    const std::string& mode,
    core::Logger& logger)
{
    if (!validateRange(threshold,
                         static_cast<int>(core::BLACK),
                         static_cast<int>(core::WHITE),
                         "threshold",
                         logger))
    {
        return false;
    }

    if (mode != "Peaks" &&
        mode != "Valleys" &&
        mode != "All")
    {
        logger.error(
            "Binary mode must be Peaks, Valleys or All");

        return false;
    }

    return true;
}

bool Validator::validateBmpWriteMode(
    const std::string& filename,
    const std::string& mode,
    core::Logger& logger)
{
    if (!validateFileName(filename, logger))
        return false;
    
    if (mode == "Full" ||
        mode == "Binary")
        return true;

    logger.error(
        "Unknown BMP mode: '" +
        mode +
        "'. Allowed values: Full, Binary");

    return false;
}

bool Validator::validateNoiseSize(
    int noisy,
    core::Logger& logger)
{
    return validateRange(noisy, 0,
                         std::numeric_limits<int>::max(),
                         "noisy", logger);
}

bool Validator::validateKMeans(
    int k,
    core::Logger& logger)
{
    return validateRange(k, 1,
                         std::numeric_limits<int>::max(),
                         "k", logger);
}

bool Validator::validateKernelSize(
    int k,
    int kernelSize,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(kernelSize,
                        1,
                        std::numeric_limits<int>::max(),
                        "kernelSize",
                        logger);

    ok &= validateKMeans(k, logger);

    return ok;
}

bool Validator::validateNavigationParameters(
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(vehicleRadius,
                        0,
                        std::numeric_limits<int>::max(),
                        "vehicleRadius",
                        logger);

    ok &= validateRange(maxSideAngle,
                        0.0,
                        90.0,
                        "maxSideAngle",
                        logger);

    ok &= validateRange(maxUpDownAngle,
                        0.0,
                        90.0,
                        "maxUpDownAngle",
                        logger);

    return ok;
}

bool Validator::validateConnectParameters(
    const DispatcherParams& params,
    const std::string& mode,
    core::Logger& logger)
{
    bool ok = true;

    ok &= validateRange(params.startPointX,
                        0,
                        params.fieldWidth - 1,
                        "startPointX",
                        logger);

    ok &= validateRange(params.startPointY,
                        0,
                        params.fieldHeight - 1,
                        "startPointY",
                        logger);

    ok &= validateRange(params.endPointX,
                        0,
                        params.fieldWidth - 1,
                        "endPointX",
                        logger);

    ok &= validateRange(params.endPointY,
                        0,
                        params.fieldHeight - 1,
                        "endPointY",
                        logger);
                        
    if (mode != "Nearest" &&
        mode != "NearestK" &&
        mode != "All")
    {
        logger.error(
            "Connect mode must be Nearest, NearestK or All");

        ok = false;
    }

    if (mode == "NearestK")
    {
        if (!validateRange(
                params.nearestVerticesCount,
                1,
                std::numeric_limits<int>::max(),
                "nearestVerticesCount",
                logger))
        {
            ok = false;
        }
    }
    
    return ok;
}

} // namespace command
