#ifndef READFILE_H
#define READFILE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "Garfield/SolidBox.hh" 


struct SolidData {
    double x, y, z; // Center position
    double dx, dy, dz; // Dimensions
};

std::vector<SolidData> loadGeometryFromFile(const std::string& filename) {

  std::vector<SolidData> solids;
  std::ifstream file(filename);
  
  if (!file.is_open()) {
    std::cout << "Error: Could not open file " << filename << std::endl;
    return solids;
  }
  
  SolidData data;
  data.dx = 0.01; // in cm
  data.dy = 0.01; // in cm
  data.dz = 0.01; // in cm
  
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    std::string token;
    
    // Read x
    std::getline(iss, token, ',');
    data.x = std::stod(token);
    
    // Read y
    std::getline(iss, token, ',');
    data.y = std::stod(token);
    
    // Read z
    std::getline(iss, token, ',');
    data.z = std::stod(token);
    
    // Read epsilon
    std::getline(iss, token, ',');
    double eps = std::stod(token);
    
    // Check if epsilon is 2.1
    if (eps == 2.1) {
      data.x /= 10; // Convert to cm
      data.y /= 10;
      data.z /= 10;
      solids.push_back(data);
    }
  }
  
  return solids;
}


// Structure to hold electric field data
struct FieldPoint {
    double x, y, z; // Position (in cm)
    double ex, ey, ez; // Electric field components (in V/cm)
};

// Function to parse the electric field file
std::vector<FieldPoint> parseFieldFile(const std::string& filename) {
    std::vector<FieldPoint> fieldData;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return fieldData;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        FieldPoint point;
        char comma; // To handle the commas in the CSV file
        if (iss >> point.x >> comma >> point.y >> comma >> point.z >> comma >>
            point.ex >> comma >> point.ey >> comma >> point.ez) {
            fieldData.push_back(point);
        }
    }

    file.close();
    return fieldData;
}

// Function to interpolate the electric field at a given point
void interpolateField(const std::vector<FieldPoint>& fieldData, double x, double y, double z,
                     double& ex, double& ey, double& ez) {
    // Simple nearest-neighbor interpolation (replace with a better method if needed)
    ex = ey = ez = 0.0;
    //double minDistance = std::numeric_limits<double>::max();
    double minDistance =0.01;
    for (const auto& point : fieldData) {
        double dx = x - point.x;
        double dy = y - point.y;
        double dz = z - point.z;
        double distance = dx * dx + dy * dy + dz * dz;
        if (distance < minDistance) {
            minDistance = distance;
            ex = point.ex;
            ey = point.ey;
            ez = point.ez;
        }
    }
}






#endif // READFILE_H
