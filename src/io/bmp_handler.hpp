#pragma once
#include <vector>
#include <string>
#include <memory>

#include "core/logger.hpp"
#include "core/constants.hpp"
#include "algorithms/gauss/pole.hpp"
#include "algorithms/gauss/gauss_builder.hpp"

namespace io {

enum class BmpWriteMode {
    Binary,
    Full
};

class BmpHandler {
private:
    core::Logger& logger;
    void logBmpHeader(const unsigned char* header) const;

public:
    BmpHandler(core::Logger& lg);

    void bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename);
    void bmp_read(algorithms::gauss::GaussBuilder& gaussBuilder, const std::string& filename,
                  std::vector<std::vector<double>>& pixelMatrix, std::unique_ptr<algorithms::gauss::Pole>& p);
};

}
