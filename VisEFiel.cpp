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
const int nx = 180, ny = 180, nz = 180;  // Grid size
const double V1 = -11000.0, V2 = -100.0;    // Plate potentials
const double epsilon_teflon = 2.1;    // Relative permittivity
//const int max_iterations = 15000;
const int max_iterations = 15000;
const double tolerance = 1e-4;

// 3D arrays for potential and permittivity
double V[nx][ny][nz]; 
double epsilon[nx][ny][nz]; 
double Ex[nx][ny][nz];
double Ey[nx][ny][nz];
double Ez[nx][ny][nz];


double grid_spacing = 0.1; // [mm/bin]

int Plate1Position = 75; // 0->Plate1Position
int TeflonLayerThick = 50; // center-TeflonLayerThick/2 -> center+TeflonLayerThick/2
int hole_radius = 13;
int hole_pitch = 50;

int ring1 = 25;
int ring2 = 50;

int ring_thz = 3;
int ting_thxy = 5;

int Plate2Position = Plate1Position+TeflonLayerThick+1; //center+TeflonLayerThick/2 

std::vector<int> planePoints = {0, nz - 1};

// Initialize grid and boundary conditions
void initialize() {

  //Centers of the system
  int cx = nx / 2, cy = ny / 2;
  
  std::cout << "initialize:: InitializeAllPotentials" << std::endl;
  ///////////////////////////////////////////////////
  //SetAllPotential to zero and permittivity to one//
  ///////////////////////////////////////////////////
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
	epsilon[i][j][k] = 1.0; 
	V[nx][ny][nz] = 0.;
      }
    }
  }

  std::cout << "initialize:: SetPlatePotentials" << std::endl;
  ///////////////////////////////////////////////////
  // Set fixed potentials///
  ///////////////////////////////////////////////////
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      V[i][j][0] = V1-3750;
      if(sqrt(pow(i - cx, 2) + pow(j - cy, 2)) > cx-5 && sqrt(pow(i - cx, 2) + pow(j - cy, 2)) < cx-5+ting_thxy){
	
	for(int k=0;k<ring_thz;k++){
	  V[i][j][ring1+k] = V1-2500;
	  V[i][j][ring2+k] = V1-1250;
	}
      }
      V[i][j][Plate1Position] = V1;
      V[i][j][Plate2Position] = V2;
      V[i][j][nz-1] = 10;
    }
  }
  
  std::cout << "initialize:: SetTeflonPermittivity" << std::endl;
  /////////////////////////////////////////////////// 
  // Define the Teflon permittivity
  /////////////////////////////////////////////////// 
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = Plate1Position; k < Plate1Position+TeflonLayerThick+1; k++) {
	epsilon[i][j][k] = epsilon_teflon;
      }
    }
  }

  /////////////////////////////////////////////////// 
  // Define the copper ring and plate permittivity
  /////////////////////////////////////////////////// 
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      epsilon[i][j][Plate1Position] = 2.5e6;
      epsilon[i][j][Plate2Position] = 2.5e6;
      epsilon[i][j][nz-1] = 2.5e6;
      epsilon[i][j][0]=2.5e6;
      for (int k=0;k<ring_thz;k++) {
	if(sqrt(pow(i - cx, 2) + pow(j - cy, 2)) > cx-5 && sqrt(pow(i - cx, 2) + pow(j - cy, 2)) < cx-5+ting_thxy){
	  epsilon[i][j][ring1+k] = 2.5e6;
	  epsilon[i][j][ring2+k] = 2.5e6;
	}
      }
    }
  }
    
  std::cout << "initialize:: SetHolesPermittivity" << std::endl;
  
  /////////////////////////////////////////////////// 
  // Define other holes in the Teflon layer (vacuum region)
  /////////////////////////////////////////////////// 
  for(int l=-1;l<=1;l++){
    for(int m=-1;m<=1;m++){
      for (int i = 0; i < nx; i++) {
	for (int j = 0; j < ny; j++) {
	  if (sqrt(pow(i - cx + l*hole_pitch, 2) + pow(j - cy+ m*hole_pitch, 2)) < hole_radius) {
	    for (int k = Plate1Position-1; k < Plate1Position+TeflonLayerThick+1+1;k++) {
	      epsilon[i][j][k] = 1.0;  // Vacuum
	    }
	  }
	}
      }
    } 
  }
  
  

}


void writePermittivityToFile(const std::string& filename) {
  // Open the file for writing
    std::ofstream outfile(filename);

    // Loop over all i, j, k
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            for (int k = 0; k < nz; ++k) {
                // Calculate the position
	      double x = (i - nx / 2) * grid_spacing / 10; //in cm
	      double y = (j - ny / 2) * grid_spacing/ 10; //in cm  
	      double z = (k - nz / 2) * grid_spacing/ 10; //in cm  
	      
	      // Write to the file
	      outfile << x << "," << y << "," << z << "\n";
            }
        }
    }

    // Close the file
    outfile.close();
    std::cout << "Data written to " << filename << std::endl;
}

void writeMeshToFile(const std::string& filename) {

  // Initialize min and max values
  double xmin = std::numeric_limits<double>::max();
  double ymin = std::numeric_limits<double>::max();
  double zmin = std::numeric_limits<double>::max();
  double xmax = std::numeric_limits<double>::lowest();
  double ymax = std::numeric_limits<double>::lowest();
  double zmax = std::numeric_limits<double>::lowest();

  // Initialize min and max values
  double xminw = std::numeric_limits<double>::max();
  double yminw = std::numeric_limits<double>::max();
  double zminw = std::numeric_limits<double>::max();
  double xmaxw = std::numeric_limits<double>::lowest();
  double ymaxw = std::numeric_limits<double>::lowest();
  double zmaxw = std::numeric_limits<double>::lowest();
  
  // Open the file for writing
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing!" << std::endl;
    return;
  }
  
  // Loop over all i, j, k
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
	// Check if epsilon[i][j][k] == 2.1
	
	// Calculate the position in cm
	double x = (i - nx / 2) * grid_spacing / 10; // in cm
	double y = (j - ny / 2) * grid_spacing / 10; // in cm
	double z = (k - nz / 2) * grid_spacing / 10; // in cm

	if (x < xminw) xminw = x;
	if (y < yminw) yminw = y;
	if (z < zminw) zminw = z;
	if (x > xmaxw) xmaxw = x;
	if (y > ymaxw) ymaxw = y;
	if (z > zmaxw) zmaxw = z;
	
	if (epsilon[i][j][k] == 2.1) {  
	  // Update min and max values
	  if (x < xmin) xmin = x;
	  if (y < ymin) ymin = y;
	  if (z < zmin) zmin = z;
	  if (x > xmax) xmax = x;
	  if (y > ymax) ymax = y;
	  if (z > zmax) zmax = z;
	  
	  // Write to the file
	  
	}
      }
    }
  }
 
  outfile << "#xminw"<< " " <<"yminw" << " " << "zminw" << " " <<"xmaxw"<< " " <<"ymaxw"<< " " <<"zmaxw"<< "\n";
  outfile << xminw<< " " <<yminw << " " << zminw << " " <<xmaxw<< " " <<ymaxw<< " " <<zmaxw<< "\n";
  outfile << "#xmin"<< " " <<"ymin" << " " << "zmin" << " " <<"xmax"<< " " <<"ymax"<< " " <<"zmax"<< "\n";
  outfile << xmin<< " " <<ymin << " " << zmin << " " <<xmax<< " " <<ymax<< " " <<zmax<< "\n";
  outfile << "#hole_rad" << " " << "pitch" << "\n"; 
  outfile << hole_radius * grid_spacing / 10 << " " << hole_pitch * grid_spacing / 10 << std::endl;
  // Close the file
  outfile.close();
  std::cout << "Data written to " << filename << std::endl;
  
  // Print min and max values
  std::cout << "xmin: " << xmin << " cm, ymin: " << ymin << " cm, zmin: " << zmin << " cm" << std::endl;
  std::cout << "xmax: " << xmax << " cm, ymax: " << ymax << " cm, zmax: " << zmax << " cm" << std::endl;
  
  
}

// Solve Laplace Equation for the potential
void solveLaplace() {
  for (int iter = 0; iter < max_iterations; iter++) {
    double max_diff = 0.0;  // Reset max difference to zero each iteration
    
    // Update interior points
    for (int i = 1; i < nx - 1; i++) {
      for (int j = 1; j < ny - 1; j++) {
	for (int k = 1; k < nz - 1; k++) {
	  // Skip fixed plate positions
	  //if ((std::find(planePoints.begin(), planePoints.end(), k) != planePoints.end() ) ||
	  //    ((k == Plate1Position || k == Plate2Position ) && sqrt(pow(i-nx/2,2) + pow(j-ny/2,2)) > hole_radius) ) continue;  

	  if(epsilon[i][j][k] == 2.5e6) continue;
	  
	  double numerator =
	    epsilon[i + 1][j][k] * V[i + 1][j][k] +
	    epsilon[i - 1][j][k] * V[i - 1][j][k] +
	    epsilon[i][j + 1][k] * V[i][j + 1][k] +
	    epsilon[i][j - 1][k] * V[i][j - 1][k] +
	    epsilon[i][j][k + 1] * V[i][j][k + 1] +
	    epsilon[i][j][k - 1] * V[i][j][k - 1];
	  
	  double denominator =
	    epsilon[i + 1][j][k] +
	    epsilon[i - 1][j][k] +
	    epsilon[i][j + 1][k] +
	    epsilon[i][j - 1][k] +
	    epsilon[i][j][k + 1] +
	    epsilon[i][j][k - 1];

	  double V_new = numerator / denominator;
	  
	  max_diff = std::max(max_diff, fabs(V_new - V[i][j][k]));
	  V[i][j][k] = V_new;
	}
      }
    }

        
    if (iter % 100 == 0) {
      std::cout << "solveLaplace::Iteration: " << iter << ", max diff: " << max_diff << std::endl;
    }
    if (max_diff < tolerance) {
      std::cout << "Converged in " << iter << " iterations.\n";
      break; 
    } 
  }
}

// Visualize the potential using ROOT
void visualize() {
    TCanvas *c1 = new TCanvas("c1", "Electric Potential", 800, 600);
    TH3F *h3 = new TH3F("h3", "Electric Potential;X;Y;Z", nx, 0, nx, ny, 0, ny, nz, 0, nz);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                h3->SetBinContent(i + 1, j + 1, k + 1, V[i][j][k]);
            }
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

    for (int i = 1; i < nx - 1; i += 1) {
      for (int j = 1; j < ny - 1; j += 1) {
	for (int k = 1; k < nz - 1; k += 1) {
	  // Compute electric field components
	  Ex[i][j][k] = -(V[i + 1][j][k] - V[i - 1][j][k]) / (2 * grid_spacing);
	  Ey[i][j][k] = -(V[i][j + 1][k] - V[i][j - 1][k]) / (2 * grid_spacing);
	  Ez[i][j][k] = -(V[i][j][k + 1] - V[i][j][k - 1]) / (2 * grid_spacing);
	  
	  double magnitude = sqrt(Ex[i][j][k] * Ex[i][j][k] + Ey[i][j][k] * Ey[i][j][k] + Ez[i][j][k] * Ez[i][j][k]);
	  
	  if(j == ny/2){
	    fmagn->SetBinContent(i,k,magnitude);
	    
	    if (magnitude > 1 && i%3==0 && k%3==0) { // Avoid plotting zero fields
	      // Project the arrow onto the xz-plane
	      double x_start = i;
	      double z_start = k;
	      double x_end = x_start + Ex[i][j][k]/magnitude;
	      double z_end = z_start + Ez[i][j][k]/magnitude;
	      
	      // Add the starting point to the TGraph
	      points->SetPoint(points->GetN(), x_start, z_start);
	      
	      // Draw the arrow using TArrow
	      TArrow *arrow = new TArrow(x_start, z_start, x_end, z_end, 0.01);
	      arrow->SetLineColor(kRed); // Set arrow color
	      arrow->SetLineWidth(1);    // Set arrow thickness
	      arrow->Draw();
	      arrow_count++;
	      
	    }//chiudo if magnitude
	  }//chiudo if y plane 
	}
      }//chiudo for i j k here around
    }
    
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

void writeElectricFieldToFile(const std::string& filename) {
  // Open the file for writing
  std::ofstream outfile(filename);
  if (!outfile.is_open()) {
    std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    return;
  }
  
  // Loop over all i, j, k (excluding boundaries)
  for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int k = 1; k < nz - 1; ++k) {

	// Calculate the position
	double x = (i - nx / 2) * grid_spacing / 10; //in cm
	double y = (j - ny / 2) * grid_spacing / 10; //in cm
	double z = (k - nz / 2) * grid_spacing / 10; //in cm
	
	// Write to the file
	outfile << x << " " << y << " " << z << " " << Ex[i][j][k]*10 << " " << Ey[i][j][k]*10 << " " << Ez[i][j][k]*10 << "\n";
      }
    }
  }
  
  // Close the file
  outfile.close();
  std::cout << "Electric field data written to " << filename << std::endl;
}

// Visualize Potential Distribution
void visualizePotential(TString name,TString name_his) {
    TCanvas *c2 = new TCanvas(name, name, 800, 600);
    TH3F *h3 = new TH3F(name_his, name_his, nx, 0, nx, ny, 0, ny, nz, 0, nz);

    std::cout <<"visualizePotentials::SettingPlot" << std::endl; 
    
    for (int i = 0; i < nx; i+=4) {
        for (int j = 0; j < ny; j+=4) {
            for (int k = 0; k < nz; k++) {
	      h3->SetBinContent(i + 1, j + 1, k + 1, V[i][j][k]);
            }
        }
    }

    h3->Draw("COLZ");
    c2->Update();
}

// Visualize Potential Distribution                                                                                                                                                                                                                       
void visualizePermittivity(TString name,TString name_his) {
    TCanvas *c4 = new TCanvas(name, name, 800, 600);
    TH3F *h4 = new TH3F(name_his, name_his, nx, 0, nx, ny, 0, ny, nz, 0, nz);

    std::cout <<"visualizePermittivity::SettingPlot" << std::endl; 
    for (int i = 0; i < nx; i+=4) {
        for (int j = 0; j < ny; j+=4) {
            for (int k = 0; k < nz; k++) {
	      if(epsilon[i][j][k]>1)  h4->SetBinContent(i + 1, j + 1, k + 1, epsilon[i][j][k]);
            }
        }
    }
    h4->SetMinimum(1);
    h4->SetMaximum(3);
    h4->Draw("COLZ");
    c4->Update();
}


int main(int argc, char** argv){

  gStyle->SetOptStat(00000);
  TApplication theApp("App",&argc, argv);
  
  initialize();
  //visualizePotential("Initial_Potential","Plot_Initial_Potential");
  //visualizePermittivity("Permittivity","Permittivity_Plot");
  //writeMeshToFile("OutFileField/TeflonMesh.txt");
  solveLaplace();
  //visualizePotential("Final_Potential","Plot_Final_Potential");
  //visualizePotential();
  visualizeElectricFieldXZProjection();
  writeElectricFieldToFile("OutFileField/EFieldFile_4mesh.txt");
  //theApp.Run();
}
