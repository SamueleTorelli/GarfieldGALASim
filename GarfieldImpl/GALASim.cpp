#include <iostream>
#include <fstream>
#include <vector>

#include "Garfield/SolidBox.hh"
#include "Garfield/SolidTube.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/GeometrySimple.hh"
#include <Garfield/ComponentFieldMap.hh>
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/ViewGeometry.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/Medium.hh"
#include <Garfield/ComponentUser.hh>
#include <Garfield/AvalancheMicroscopic.hh>
#include <Garfield/ComponentGrid.hh>
#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"

#include "TApplication.h"
#include "TPolyLine3D.h"
#include <TH2D.h>

#include "Readfile.h"
#include "AnalyseAvalanche.h"

using namespace Garfield;

std::ofstream outFile; // Global file pointer

int ev;

void openFile() {
    outFile.open("example.txt", std::ios::out | std::ios::app); // Open for writing, append mode
    if (!outFile) {
        std::cerr << "Error opening file!" << std::endl;
    }
}

void writeToFile(int ev,double x, double y, double z, double t) {

  outFile << ev << " " << x << " " << y << " " << z << " " << t << "\n" ;
   
}

void closeFile() {
    if (outFile.is_open()) {
        outFile.close();
    }
}


void userHandleStep(double x, double y, double z, double t,double e, double dx, double dy, double dz, bool hole){
  if(z>0.4)    std::cout << x << "\t" << y << "\t" << z << "\t" << e << "\t" << dx << "\t" << dy << "\t" << dz << "\t" << hole << "\n";
}

void userTracking(double x, double y, double z, double t,
                int type, int level, Garfield::Medium* /*1e-2*/) {
    // Skip inelastic collisions that are not excitations.
    if (type != 4) return; //save tracks ONLY of	inelastic collisions that have produced	excitations.
    
    std::cout<< ev << " " << x << " " << y << " " << z << " " << t << std::endl;
    writeToFile(ev,x,y,z,t);
}



int main() {

  TApplication app("app", nullptr, nullptr);
    
  MediumMagboltz gas;
  gas.LoadGasFile("/Users/samuele/Documents/Postdoc/GALA/GarfieldGALADSim/GasFile/gases/xe_5bar.gas");
  gas.SetTemperature(293.15); // 20°C
  gas.SetPressure(760.*5);      // 1 atm
  gas.Initialise(true); // Enable Penning transfer
  gas.EnableDrift(true);
  
  // Define a Teflon medium
  Medium teflon;
  teflon.SetDielectricConstant(2.1);
  teflon.EnableDrift(false);
  
  // Create a geometry
  GeometrySimple geometry;
  
  // Variables to store values from the file
  double xminw, yminw, zminw, xmaxw, ymaxw, zmaxw;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  double hole_rad, pitch;
  
  std::ifstream infile("/Users/samuele/Documents/Postdoc/GALA/EFSimulation/OutFileField/TeflonMesh.txt"); // Change this to your filename
  if (!infile.is_open()) {
    std::cerr << "Error: Could not open file!" << std::endl;
    return 1;
  }
  
  std::string line;
  while (std::getline(infile, line)) {
    if (line.empty() || line[0] == '#') continue; // Skip comment lines
    
    std::istringstream iss(line);
    if (!(xminw || xmin)) {  // First non-comment line → xminw values
      iss >> xminw >> yminw >> zminw >> xmaxw >> ymaxw >> zmaxw;
    } else if (!(xmin)) {  // Second non-comment line → xmin values
      iss >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax;
    } else {  // Third non-comment line → hole_rad and pitch
      iss >> hole_rad >> pitch;
    }
  }

  zmin = zmin;
  zmax = zmax;
  
  // Close the file
  infile.close();
  std::cout << "Gas box: (" << xminw << ", " << yminw << ", " << zminw << ") to ("
	    << xmaxw << ", " << ymaxw << ", " << zmaxw << ")\n";
  std::cout << "Teflon box: (" << xmin << ", " << ymin << ", " << zmin << ") to ("
	    << xmax << ", " << ymax << ", " << zmax << ")\n";
  std::cout << "Hole radius: " << hole_rad << ", Pitch: " << pitch << "\n";

      // Create the Teflon box
  double teflonCenterX = (xmin + xmax) / 2;
  double teflonCenterY = (ymin + ymax) / 2;
  double teflonCenterZ = (zmin + zmax) / 2;
  double teflonSizeX = (xmax - xmin) / 2;
  double teflonSizeY = (ymax - ymin) / 2;
  double teflonSizeZ = (zmax - zmin) / 2;
  
  // Create 9 cylindrical holes in the Teflon box
  std::vector<std::pair<double, double>> holePositions = {
    {0.0, 0.0},          // Center
    {-pitch, 0.0},       // Left
    {pitch, 0.0},        // Right
    {0.0, pitch},        // Top
    {0.0, -pitch},       // Bottom
    {-pitch, pitch},     // Top-left
    {pitch, pitch},     // Top-right
    {-pitch, -pitch},   // Bottom-left
    {pitch, -pitch}     // Bottom-right
  };

  for (const auto& pos : holePositions) {
    double x = pos.first;
    double y = pos.second;
    double z = teflonCenterZ; // Center of the hole in z
    double height = zmax - zmin; // Height of the hole

    auto hole = new Garfield::SolidTube(x, y, z, hole_rad, height/2+0.02);
    geometry.AddSolid(hole, &gas); // Assign gas medium to the hole
  }

  auto teflonBox = new Garfield::SolidBox(teflonCenterX, teflonCenterY, teflonCenterZ, teflonSizeX, teflonSizeY, teflonSizeZ+0.01);
  geometry.AddSolid(teflonBox, &teflon);
  
      // Create the gas box
  double gasCenterX = (xminw + xmaxw) / 2;
  double gasCenterY = (yminw + ymaxw) / 2;
  double gasCenterZ = (zminw + zmaxw) / 2;
  double gasSizeX = (xmaxw - xminw) / 2;
  double gasSizeY = (ymaxw - yminw) / 2;
  double gasSizeZ = (zmaxw - zminw) / 2;
  
  auto gasBox = new Garfield::SolidBox(gasCenterX, gasCenterY, gasCenterZ, gasSizeX, gasSizeY, gasSizeZ);
  geometry.AddSolid(gasBox, &gas);
  
  // Print success message
  std::cout << "Geometry created successfully!" << std::endl;
  
  // Create the field component
  ComponentGrid* field = new ComponentGrid();
  field->LoadElectricField("/Users/samuele/Documents/Postdoc/GALA/EFSimulation/OutFileField/EFieldFile_4mesh.txt", "XYZ", false, false);
  field->SetGeometry(&geometry);
  field->EnableCustomMediumReturn();
  field->Print();
  field->SetMedium(&gas);

  // Create a sensor and associate the field component
  Sensor sensor;
  sensor.AddComponent(field);
  sensor.EnableDebugging();
  // Set the sensor region to the amplification volume
  sensor.SetArea(-0.4, -0.4, -0.3,
		  0.4, 0.4, 0.5);

  new TCanvas();
  TH2D* h = new TH2D("EF","EF",int((xmaxw-xminw)/0.01),xminw,xmaxw,int((zmaxw-zminw)/0.01),zminw,zmaxw);
  
  for(double x =  xminw+0.01; x <xmaxw-0.01; x += 0.001){
    for (double z = zminw; z <zmaxw; z += 0.001) {
      double ex, ey, ez;
      
      ex = field->ElectricField(x, 0.0, z)[0];
      ey = field->ElectricField(x, 0.0, z)[1];
      ez = field->ElectricField(x, 0.0, z)[2];
      
      double magn = sqrt(ex*ex+ey*ey+ez*ez);
      //if(z>0.2 && z<0.4 && magn < 1) std::cout <<"Magn 0 at (x,z) = " << x << " " << z << std::endl;
      h->SetBinContent(h->GetXaxis()->FindBin(x),h->GetYaxis()->FindBin(z),magn);
    }
  }

  for (double z = zminw; z <zmaxw; z += 0.001) {
    Medium* m = sensor.GetMedium(-0.25, 0, z);
    if (!m) {
      std::cout << "Warning: Electron is in an undefined region!" << std::endl;
    } else {
      std::cout << "Warning: Medium driftable? " << m->GetName() << " " << m->IsDriftable() << std::endl;
    } 
  }
  
  h->Draw("COLZ");
  
  // Create an avalanche microscopic object for electron tracking
  AvalancheMicroscopic avalanche;
  avalanche.SetSensor(&sensor); // Associate the avalanche with the sensor
  //Enable debugging output
  //avalanche.EnableDebugging();
  avalanche.EnableRKNSteps(false);

  avalanche.SetRKNTolerance(1e-5,1e-2);
  //avalanche.SetUserHandleStep(userHandleStep);
  // Set up visualization of the electron drift

  ViewDrift viewDrift;
  avalanche.EnablePlotting(&viewDrift);
  avalanche.EnableExcitationMarkers(true);
  avalanche.SetElectronTransportCut(0.00001);
  avalanche.SetUserHandleInelastic(userTracking);

  bool debug = false;
  
  if(!debug){
    openFile();
    
    std::cout << "Drifting medium: "<< gas.IsDriftable() << " " << teflon.IsDriftable() << "\n"; 
    
    for(int i=0; i<10; i++){
      // Set the initial position of the electron (in cm)
      double x0 = rand()/(double)RAND_MAX*0.4-0.2;
      double y0 = rand()/(double)RAND_MAX*0.4-0.2;
      double z0 = -0.25; // Start at the bottom of the volume
      std::cout << i <<" e- starting at " << x0 << " " << y0 << " " << z0 << std::endl;
      double t0 = 0.0; // Initial time (in ns)
      double e0 = 0.01; // Initial energy of the electron (in eV)
      
      // Simulate the electron drift within the amplification volume
      avalanche.DriftElectron(x0, y0, z0, t0, e0);
      ev=i;
      AnalyseAval(avalanche,field);
      
    }
  
    closeFile();
  }
  viewDrift.Plot();

  // Visualize the geometry
  ViewGeometry view;
  view.SetGeometry(&geometry);
  view.Plot();
    
  app.Run();
  return 0;
}
