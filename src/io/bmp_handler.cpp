#include "io/bmp_handler.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdint>

namespace io {

BmpHandler::BmpHandler(core::Logger& lg) : logger(lg) {}

    void BmpHandler::logBmpHeader(const unsigned char* header) const {
        std::ostringstream oss;
        oss << "[BmpHandler] BMP Header Info:\n"
            << "  File Size: " << *reinterpret_cast<const uint32_t*>(&header[2]) << " bytes\n"
            << "  Image Size: " << *reinterpret_cast<const uint32_t*>(&header[34]) << " bytes\n"
            << "  Dimensions: " << *reinterpret_cast<const int32_t*>(&header[18]) 
            << "x" << *reinterpret_cast<const int32_t*>(&header[22]) << "\n"
            << "  Bits Per Pixel: " << *reinterpret_cast<const uint16_t*>(&header[28]);
        logger.debug(oss.str());
    }

    void BmpHandler::bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename) {
        logger.trace("[BmpHandler::bmp_write] Starting BMP write operation");
        
        if (pixelMatrix.empty() || pixelMatrix[0].empty()) {
            logger.error("[BmpHandler::bmp_write] Empty pixel matrix provided");
            return;
        }

        int width = pixelMatrix[0].size();
        int height = pixelMatrix.size();
        int padding = (4 - (width * 3) % 4) % 4;

        logger.debug(std::string("[BmpHandler::bmp_write] Creating BMP file: ") + filename + 
                   " with dimensions " + std::to_string(width) + "x" + std::to_string(height) +
                   ", padding: " + std::to_string(padding) + " bytes");

        std::ofstream bmpFile(filename, std::ios::binary);
        if (!bmpFile) {
            logger.error(std::string("[BmpHandler::bmp_write] Failed to create BMP file: ") + filename);
            return;
        }

        // Prepare and write header
    unsigned char bmpHeader[54] = {
        'B', 'M', // Identifier
        0, 0, 0, 0, // Size of file (will be set later)
        0, 0, 0, 0, // Reserved
        54, 0, 0, 0, // Header size
        40, 0, 0, 0, // Info header size
        0, 0, 0, 0, // Width (will be set later)
        0, 0, 0, 0, // Height (will be set later)
        1, 0, // Number of color planes
        24, 0, // Bits per pixel
        0, 0, 0, 0, // Compression
        0, 0, 0, 0, // Image size (will be set later)
        0x13, 0x0B, 0, 0, // Horizontal resolution
        0x13, 0x0B, 0, 0, // Vertical resolution
        0, 0, 0, 0, // Number of colors in palette
        0, 0, 0, 0  // Important colors
    };
    // Set width and height in header
    bmpHeader[18] = (width & 0xFF);
    bmpHeader[19] = (width >> 8) & 0xFF;
    bmpHeader[20] = (width >> 16) & 0xFF;
    bmpHeader[21] = (width >> 24) & 0xFF;
    bmpHeader[22] = (height & 0xFF);
    bmpHeader[23] = (height >> 8) & 0xFF;
    bmpHeader[24] = (height >> 16) & 0xFF;
    bmpHeader[25] = (height >> 24) & 0xFF;
    // Write header
    logger.trace("[BmpHandler::bmp_write] Writing BMP header");
        bmpFile.write(reinterpret_cast<char*>(bmpHeader), 54);

        // Write pixel data
        size_t pixels_written = 0;
        size_t clamped_pixels = 0;
        logger.trace("[BmpHandler::bmp_write] Writing pixel data (bottom-to-top)");
        
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                double pixelValue = pixelMatrix[y][x];
                if (pixelValue < core::BLACK || pixelValue > core::WHITE) {
                    clamped_pixels++;
                    pixelValue = std::clamp(pixelValue,
                       static_cast<double>(core::BLACK),
                       static_cast<double>(core::WHITE));
                }
                unsigned char color = static_cast<unsigned char>(pixelValue);
                bmpFile.put(color).put(color).put(color);
                pixels_written++;
            }
            // Write padding
            for (int p = 0; p < padding; ++p) {
                bmpFile.put(0);
            }
        }

        bmpFile.close();
        logger.info(std::string("[BmpHandler::bmp_write] BMP file successfully written\n") +
                  "  Total pixels: " + std::to_string(pixels_written) + "\n" +
                  "  Clamped pixels: " + std::to_string(clamped_pixels) + "\n" +
                  "  File size: ~" + std::to_string(54 + (width*3 + padding)*height) + " bytes");
    }
    
   void BmpHandler::bmp_read(algorithms::gauss::GaussBuilder& gaussBuilder, const std::string& filename, 
                 std::vector<std::vector<double>>& pixelMatrix, std::unique_ptr<algorithms::gauss::Pole>& p) {
        logger.trace(std::string("[BmpHandler::bmp_read] Starting BMP read operation: ") + filename);
        
        std::ifstream bmpFile(filename, std::ios::binary);
        if (!bmpFile) {
            logger.error(std::string("[BmpHandler::bmp_read] Failed to open BMP file: ") + filename);
            return;
        }

        // Read header
        unsigned char header[54];
        bmpFile.read(reinterpret_cast<char*>(header), 54);
        logBmpHeader(header);

        int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
        int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
        int padding = (4 - (width * 3) % 4) % 4;

        logger.debug(std::string("[BmpHandler::bmp_read] Image dimensions: ") + 
                   std::to_string(width) + "x" + std::to_string(height) +
                   ", padding: " + std::to_string(padding) + " bytes");

        // Initialize field
        logger.trace("[BmpHandler::bmp_read] Initializing field through GaussBuilder");
        gaussBuilder.init(height, width, p);
        pixelMatrix.resize(height, std::vector<double>(width));

        // Read pixel data
        size_t pixels_read = 0;
        double min_val = core::WHITE, max_val = core::BLACK;
        logger.trace("[BmpHandler::bmp_read] Reading pixel data (bottom-to-top)");
        
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                unsigned char color = bmpFile.get();
                bmpFile.get(); // G
                bmpFile.get(); // R
                
                double value = static_cast<double>(color);
                pixelMatrix[y][x] = value;
                p->field[y][x] = value;
                
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
                pixels_read++;
            }
            bmpFile.ignore(padding);
        }

        bmpFile.close();
        logger.info(std::string("[BmpHandler::bmp_read] BMP file successfully loaded\n") +
                  "  Pixels read: " + std::to_string(pixels_read) + "\n" +
                  "  Value range: [" + std::to_string(min_val) + ".." + 
                  std::to_string(max_val) + "]");
    }
}
