#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include <TCanvas.h>
#include <TH3F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TArrow.h>
#include <TStyle.h>
#include <TApplication.h>

using namespace std;

// Simulation parameters
const int nx = 300, nz = 1000;  // Grid size
const double V1 = -11000.0, V2 = -100.0;    // Plate potentials
const double epsilon_teflon = 2.1;    // Relative permittivity
const double epsilon_copper = 2e5;

const int max_iterations = 150000;
const double tolerance = 0.0001;

// 3D arrays for potential and permittivity
// Define 3D vectors for V and epsilon
std::vector<std::vector<double>> V(nx,std::vector<double>(nz, 0.0));
std::vector<std::vector<double>> epsilon(nx,std::vector<double>(nz, 1.0)); // Default value for epsilon (e.g., 1.0)

double grid_spacing = 0.025; // [mm/bin]

int Plate1Position = 400; // 0->Plate1Position
int TeflonLayerThick = 200; // center-TeflonLayerThick/2 -> center+TeflonLayerThick/2
int hole_radius = 100/2;
int hole_pitch = 200;

int Plate2Position = Plate1Position+TeflonLayerThick+1; //center+TeflonLayerThick/2 

// Initialize grid and boundary conditions

// Initialize grid and boundary conditions
void initialize() {

  //Centers of the system
  int cx = nx / 2;
  
  std::cout << "initialize:: SetPlatePotentials" << std::endl;
  ///////////////////////////////////////////////////
  // Set fixed potentials///
  ///////////////////////////////////////////////////
  for (int i = 0; i < nx; i++) {
    V[i][1] = V1-500;
    V[i][Plate1Position] = V1;
    V[i][Plate2Position] = V2;
    V[i][nz-1] = 50;
  }
  
  std::cout << "initialize:: SetTeflonPermittivity" << std::endl;
  /////////////////////////////////////////////////// 
  // Define the Teflon permittivity
  /////////////////////////////////////////////////// 
  for (int i = 0; i < nx; i++) {
      for (int k = Plate1Position; k < Plate1Position+TeflonLayerThick+1; k++) {
	epsilon[i][k] = epsilon_teflon;
      }
  }

  /////////////////////////////////////////////////// 
  // Define the copper ring and plate permittivity
  /////////////////////////////////////////////////// 
  for (int i = 0; i < nx; i++) {
    epsilon[i][Plate1Position] = epsilon_copper;
    epsilon[i][Plate2Position] = epsilon_copper;
    epsilon[i][nz-1] = epsilon_copper;
    epsilon[i][1]=epsilon_copper;
  }
  
  std::cout << "initialize:: SetHolesPermittivity" << std::endl;
  
  /////////////////////////////////////////////////// 
  // Define other holes in the Teflon layer (vacuum region)
  /////////////////////////////////////////////////// 
  for(int l=0;l<1;l++){
    for(int m=0;m<1;m++){
      for (int i = 0; i < nx; i++) {
	if (abs(i - cx + l*hole_pitch) < hole_radius) {
	  for (int k = Plate1Position-1; k < Plate1Position+TeflonLayerThick+1+1;k++) {
	    epsilon[i][k] = 1.0;  // Vacuum
	  }
	}
      }
    } 
  }

}

void solveLaplaceSOR()
{
    // Over-relaxation parameter, typically 1 < omega < 2. Adjust for optimal performance.
    double omega = 1.7;
    int iter = 0;
    double max_diff;

    do {
      max_diff = 0.0;

      for (int k = 0; k < nz; ++k) {
	if (std::abs(epsilon[0][k] - epsilon_copper) < 1e-6)  continue;
	V[0][k]     = V[1][k];
	V[nx-1][k]  = V[nx-2][k];
      }
      
      // Loop over interior points. Adjust if your boundaries need special handling.
      for (int k = 1; k < nz-1; ++k) {   // Move k-loop to outermost
	for (int i = 1; i < nx-1; ++i) {
	  // Skip fixed copper points.
	  if (std::abs(epsilon[i][k] - epsilon_copper) < 1e-6)  continue;
          
	  // Calculate weighted contributions from the six neighbors.
	  // The weighting uses the epsilon of each neighbor.
	  
	  double weightSum = 
	    epsilon[i+1][k] + epsilon[i-1][k] +
	    epsilon[i][k+1] + epsilon[i][k-1];
	  
	  double neighborSum = 
	    epsilon[i+1][k] * V[i+1][k] +
	    epsilon[i-1][k] * V[i-1][k] +
	    epsilon[i][k+1] * V[i][k+1] +
	    epsilon[i][k-1] * V[i][k-1];
	  
	  // Compute the new potential value from the weighted average.
	  double newV = neighborSum / weightSum;
	  double diff = newV - V[i][k];
          
	  // Update using SOR.
	  V[i][k] += omega * diff;
          
	  // Track the maximum change over the grid.
	  max_diff = std::max(max_diff, std::fabs(omega * diff));
	}
      }

      
      
      
      iter++;
      // Optionally, output iteration status:
      if(iter%50==0)std::cout << "Iteration " << iter << ", max difference = " << max_diff << std::endl;
    } while (max_diff > tolerance && iter < max_iterations);
    
    std::cout << "SOR converged in " << iter << " iterations with max diff " << max_diff << std::endl;
}


// Visualize the potential using ROOT
void visualize() {
    TCanvas *c1 = new TCanvas("c1", "Electric Potential", 800, 600);
    TH2F *h3 = new TH2F("h3", "Electric Potential;X;Y;Z", nx, 0, nx, nz, 0, nz);

    for (int i = 0; i < nx; i++) {
      for (int k = 0; k < nz; k++) {
	h3->SetBinContent(i , k , V[i][k]);
      }
    }

    h3->Draw("COLZ");
    c1->Update();
}


void visualizeElectricFieldXZProjection() {
    TCanvas *c1 = new TCanvas("c1", "Electric Field XZ Projection", 800, 600);
    c1->SetRightMargin(0.15); // Adjust margin for better visualization
          
    // Create a TH2F to define the axes and frame
    TH2F *frame = new TH2F("frame", "Electric Field XZ Projection; X; Z", 
                            nx, 0, nx, nz, 0, nz); // Adjust binning and range as needed
    frame->SetStats(0); // Disable statistics box
    frame->Draw(); // Draw the frame to set up the axes

    TH2F* fmagn = new TH2F("fmagn", "EFMagnitude",nx, 0, nx, nz, 0, nz);
    fmagn->Draw("COLZ");

    // Create a TGraph to plot the starting points
    TGraph *points = new TGraph();
    points->SetMarkerStyle(20); // Set marker style (20 = small circle)
    points->SetMarkerSize(0.2); // Set marker size
    points->SetMarkerColor(kBlue); // Set marker color

    int arrow_count = 0;
    std::cout << "visualizeElectricFieldXZProjection::Calculating electric fields" << std::endl;

    double Ex,Ez;
    
    for (int i = 2; i < nx - 2; i += 1) {
      for (int k = 2; k < nz - 2; k += 1) {
	
	// Compute electric field components using fourth-order differences
	Ex = -(-V[i+2][k] + 8.0 * V[i+1][k] - 8.0 * V[i-1][k] + V[i-2][k]) / (12.0 * grid_spacing);
	Ez = -(-V[i][k+2] + 8.0 * V[i][k+1] - 8.0 * V[i][k-1] + V[i][k-2]) / (12.0 * grid_spacing);
	
	double magnitude = sqrt(Ex * Ex +  Ez * Ez);
	
	fmagn->SetBinContent(i,k,magnitude);
	
	if (magnitude > 1 && i%5==0 && k%5==0) { // Avoid plotting zero fields
	  // Project the arrow onto the xz-plane
	  double x_start = i;
	  double z_start = k;
	  double x_end = x_start + Ex/magnitude;
	  double z_end = z_start + Ez/magnitude;
	  
	  // Add the starting point to the TGraph
	  points->SetPoint(points->GetN(), x_start, z_start);
	  
	  // Draw the arrow using TArrow
	  TArrow *arrow = new TArrow(x_start, z_start, x_end, z_end, 0.01);
	  arrow->SetLineColor(kRed); // Set arrow color
	  arrow->SetLineWidth(1);    // Set arrow thickness
	  arrow->Draw();
	  arrow_count++;
	  
	}//chiudo if magnitude
      }
    }//chiudo for i j k here around

    
    // Draw the points on the canvas
    //points->Draw("P SAME"); // "P" for points, "SAME" to overlay on existing canvas
    
    std::cout << "Drew " << arrow_count << " arrows in the XZ plane." << std::endl;
    c1->Update();
    /*
    TCanvas *c2 = new TCanvas("c7", "Electric Field magn XZ Projection", 800, 600);
    fmagn->Draw("COLZ");
    c2->Update();
    */
}


// Visualize Potential Distribution
void visualizePotential(TString name,TString name_his) {
    TCanvas *c2 = new TCanvas(name, name, 800, 600);
    TH2F *h3 = new TH2F(name_his, name_his, nx, 0, nx, nz, 0, nz);

    std::cout <<"visualizePotentials::SettingPlot" << std::endl; 
    
    for (int i = 0; i < nx; i+=1) {
            for (int k = 0; k < nz; k++) {
	      h3->SetBinContent(i, k , V[i][k]);
            }
    }

    h3->Draw("COLZ");
    c2->Update();
}

// Visualize Potential Distribution                                                                                                                                                                                                                       
void visualizePermittivity(TString name,TString name_his) {
    TCanvas *c4 = new TCanvas(name, name, 800, 600);
    TH2F *h4 = new TH2F(name_his, name_his, nx, 0, nx, nz, 0, nz);

    std::cout <<"visualizePermittivity::SettingPlot" << std::endl; 
    for (int i = 0; i < nx; i+=1) {
      for (int k = 0; k < nz; k++) {
	if(epsilon[i][k]>1)  h4->SetBinContent(i, k, epsilon[i][k]);
      }
    }

    h4->SetMinimum(1);
    h4->SetMaximum(3);
    h4->Draw("COLZ");
    c4->Update();
}

void writeMeshToFile(const std::string& filename) {

  // Initialize min and max values
  double xmin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();
  double zmax = std::numeric_limits<double>::lowest();

  // Initialize min and max values
  double xminw = std::numeric_limits<double>::max();
  double zminw = std::numeric_limits<double>::max();
  double xmaxw = std::numeric_limits<double>::lowest();
  double zmaxw = std::numeric_limits<double>::lowest();
  
  // Open the file for writing
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
    return;
  }
  
  // Loop over all i, j, k
  for (int i = 0; i < nx; ++i) {
    for (int k = 0; k < nz; ++k) {
      // Check if epsilon[i][j][k] == 2.1
      
      // Calculate the position in cm
      double x = (i - nx / 2) * grid_spacing / 10; // in cm
      double z = (k - nz / 2) * grid_spacing / 10; // in cm
      
      if (x < xminw) xminw = x;
      if (z < zminw) zminw = z;
      if (x > xmaxw) xmaxw = x;
      if (z > zmaxw) zmaxw = z;
      
      if (epsilon[i][k] == 2.1) {  
	// Update min and max values
	if (x < xmin) xmin = x;
	if (z < zmin) zmin = z;
	if (x > xmax) xmax = x;
	if (z > zmax) zmax = z;
	
	// Write to the file
	
      }
    }
  }
  
  outfile << "#xminw"<< " " <<"yminw" << " " << "zminw" << " " <<"xmaxw"<< " " <<"ymaxw"<< " " <<"zmaxw"<< "\n";
  outfile << xminw<< " " <<xminw << " " << zminw << " " <<xmaxw<< " " <<xmaxw<< " " <<zmaxw<< "\n";
  outfile << "#xmin"<< " " <<"ymin" << " " << "zmin" << " " <<"xmax"<< " " <<"ymax"<< " " <<"zmax"<< "\n";
  outfile << xmin<< " " <<xmin << " " << zmin << " " <<xmax<< " " <<xmax<< " " <<zmax<< "\n";
  outfile << "#hole_rad" << " " << "pitch" << "\n"; 
  outfile << hole_radius * grid_spacing / 10 << " " << hole_pitch * grid_spacing / 10 << std::endl;
  // Close the file
  outfile.close();
  std::cout << "Data written to " << filename << std::endl;
  
  // Print min and max values
  std::cout << "xmin: " << xmin << " cm, ymin: " << xmin << " cm, zmin: " << zmin << " cm" << std::endl;
  std::cout << "xmax: " << xmax << " cm, ymax: " << xmax << " cm, zmax: " << zmax << " cm" << std::endl;
  
  
}


void writeElectricFieldToFile2D(const std::string& filename) {
  // Open the file for writing
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    return;
  }

  double Ex, Ez;
  
  // Loop over all i, j, k (excluding boundaries)
  for (int i = 1; i < nx - 1; ++i) {
    for (int k = 1; k < nz - 1; ++k) {	
      
      Ex = -(-V[i+2][k] + 8.0 * V[i+1][k] - 8.0 * V[i-1][k] + V[i-2][k]) / (12.0 * grid_spacing);
      Ez = -(-V[i][k+2] + 8.0 * V[i][k+1] - 8.0 * V[i][k-1] + V[i][k-2]) / (12.0 * grid_spacing);
	
	// Calculate the position
      double x = (i - nx / 2) * grid_spacing / 10; //in cm
      double z = (k - nz / 2) * grid_spacing / 10; //in cm
      
      // Write to the file
      outfile << x <<  " " << z << " " << Ex*10  << " " << Ez*10 << "\n";
    }
  }
  

  // Close the file
  outfile.close();
  std::cout << "Electric field data written to " << filename << std::endl;
}


void writeElectricField3D(const std::string& filename, int n_phi = 72) {
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    return;
  }

  double Ex_r, Ez, Er;

  // Azimuthal step
  double dphi = 2.0 * M_PI / n_phi;

  for (int i = 2; i < nx - 2; ++i) {
    for (int k = 2; k < nz - 2; ++k) {

      // Radial field in x-z plane
      Er = -(-V[i+2][k] + 8.0 * V[i+1][k] - 8.0 * V[i-1][k] + V[i-2][k]) / (12.0 * grid_spacing);
      Ez = -(-V[i][k+2] + 8.0 * V[i][k+1] - 8.0 * V[i][k-1] + V[i][k-2]) / (12.0 * grid_spacing);

      double r = (i - nx / 2) * grid_spacing / 10; // in cm
      double z = (k - nz / 2) * grid_spacing / 10; // in cm

      // Loop over azimuthal angle phi
      for (int p = 0; p < n_phi; ++p) {
        double phi = p * dphi;

        // Cylindrical to Cartesian
        double x = r * cos(phi);
        double y = r * sin(phi);

        double Ex = Er * cos(phi);  // E_r in x
        double Ey = Er * sin(phi);  // E_r in y

        outfile << x << " " << y << " " << z << " " << Ex * 10 << " " << Ey * 10 << " " << Ez * 10 << "\n";
      }
    }
  }

  outfile.close();
  std::cout << "3D Electric field data with cylindrical symmetry written to " << filename << std::endl;
}

void writeElectricField3DXYZ(const std::string& filename) {
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    return;
  }

  std::cout << "writeElectricField3DXYZ::Writing EF " << std::endl;
  
  int ny=nx;
  int n_xy = nx;
  
  // Determine simulation radial range.
  // We use i=2 to i=nx-2, so the maximum radial distance (in cm) is:
  double r_max = ((nx/2) - 2) * grid_spacing / 10.0;

  // Define the output (x,y) grid.
  // Here we assume that the simulation is centered at x = 0 and y = 0.
  // We choose x and y to span at least [-r_max, r_max].
  double x_min = -r_max, x_max = r_max;
  double y_min = -r_max, y_max = r_max;

  // Loop over the output Cartesian grid in x and y, and then over z.
  // The z range is taken from the simulation points (we use i from 2 to nx-2 for r and k from 2 to nz-2 for z).
  for (int ix = 0; ix < n_xy; ++ix) {
    double x = x_min + (x_max - x_min) * ix / (n_xy - 1);
    for (int iy = 0; iy < n_xy; ++iy) {
      double y = y_min + (y_max - y_min) * iy / (n_xy - 1);

      // Determine the corresponding radial coordinate.
      double r = std::sqrt(x*x + y*y);
      // Determine the azimuthal angle (for converting from the radial field to Cartesian components).
      double phi = std::atan2(y, x);

      // Loop over the z-axis (using the same sampling as in the simulation).
      for (int k = 2; k < nz - 2; ++k) {
        double z = (k - nz / 2) * grid_spacing / 10.0; // in cm

        double Ex, Ey, Ez;
        if (r > r_max) {
          // Outside the available (cylindrical) simulation region.
          Ex = 0.0;
          Ey = 0.0;
          Ez = 0.0;
        } else {
          // Map the physical r back to a simulation index.
          // In the simulation, r_sim = (i - nx/2)*grid_spacing/10, so:
          double i_real = r * 10.0 / grid_spacing + nx/2;
          // Ensure we are within the valid range:
          i_real = std::clamp(i_real, 2.0, double(nx - 3));
          
          int i_low = int(std::floor(i_real));
          int i_high = int(std::ceil(i_real));

          double r_low = (i_low - nx/2) * grid_spacing / 10.0;
          double r_high = (i_high - nx/2) * grid_spacing / 10.0;

          // Compute the radial derivative (electric field magnitude in the radial direction)
          // using the finite-difference formula at the two indices.
          double Er_low = -(-V[i_low+2][k] + 8.0 * V[i_low+1][k] - 8.0 * V[i_low-1][k] + V[i_low-2][k])
                          / (12.0 * grid_spacing);
          double Er_high = -(-V[i_high+2][k] + 8.0 * V[i_high+1][k] - 8.0 * V[i_high-1][k] + V[i_high-2][k])
                           / (12.0 * grid_spacing);

          // Linear interpolation in the radial direction
          double Er_interp = (r_high != r_low)
            ? Er_low + (Er_high - Er_low) * (r - r_low) / (r_high - r_low)
            : Er_low;

          // For E_z, use the finite difference along z at the lower index.
          double Ez_val = -(-V[i_low][k+2] + 8.0 * V[i_low][k+1] - 8.0 * V[i_low][k-1] + V[i_low][k-2])
                          / (12.0 * grid_spacing);

          // Convert the radial field Er to Cartesian components.
          Ex = Er_interp * std::cos(phi);
          Ey = Er_interp * std::sin(phi);
          Ez = Ez_val;
        }

        // The simulation multiplies field values by 10 (to convert to V/cm) as in your original code.
        outfile << x << " " << y << " " << z << " "
                << Ex * 10 << " " << Ey * 10 << " " << Ez * 10 << "\n";
      } // end z-loop
    } // end y-loop
  } // end x-loop

  outfile.close();
  std::cout << "3D Electric field data written to " << filename
            << " with " << n_xy << " bins in x and y and 3D volume filled (field=0 outside cylinder)" << std::endl;
}




int main(int argc, char** argv){

  gStyle->SetOptStat(00000);
  TApplication theApp("App",&argc, argv);
  
  initialize();
  visualizePotential("Initial_Potential","Plot_Initial_Potential");
  visualizePermittivity("Permittivity","Permittivity_Plot");
  //  writeMeshToFile("OutFileField/TeflonMesh.txt");

  solveLaplaceSOR();
  
  visualizePotential("Final_Potential","Plot_Final_Potential");
  
  visualizeElectricFieldXZProjection();
  //writeElectricField3D("OutFileField/field_test.txt");
  //  writeElectricField3DXYZ("OutFileField/field_testXYZ.txt");
  theApp.Run();
}
