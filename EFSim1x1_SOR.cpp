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
const int nx = 300, ny = 300, nz = 1000;  // Grid size
const double V1 = -11000.0, V2 = -100.0;    // Plate potentials
const double epsilon_teflon = 2.1;    // Relative permittivity
const double epsilon_copper = 10;
//const int max_iterations = 15000;
const int max_iterations = 15000;
const double tolerance = 5;

// 3D arrays for potential and permittivity
// Define 3D vectors for V and epsilon
std::vector<std::vector<std::vector<double>>> V(nx, 
    std::vector<std::vector<double>>(ny, 
    std::vector<double>(nz, 0.0)));

std::vector<std::vector<std::vector<double>>> epsilon(nx, 
    std::vector<std::vector<double>>(ny, 
    std::vector<double>(nz, 1.0))); // Default value for epsilon (e.g., 1.0)

double grid_spacing = 0.025; // [mm/bin]

int Plate1Position = 400; // 0->Plate1Position
int TeflonLayerThick = 200; // center-TeflonLayerThick/2 -> center+TeflonLayerThick/2
int hole_radius = 100/2;
int hole_pitch = 200;

int Plate2Position = Plate1Position+TeflonLayerThick+1; //center+TeflonLayerThick/2 

// Initialize grid and boundary conditions
void initialize() {

  //Centers of the system
  int cx = nx / 2, cy = ny / 2;
  
  std::cout << "initialize:: InitializeAllPotentials" << std::endl;
  ///////////////////////////////////////////////////
  //SetAllPotential to zero and permittivity to one//
  ///////////////////////////////////////////////////
  /*for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      for (int k = 0; k < nz; ++k) {
	epsilon[i][j][k] = 1.0; 
	V[i][j][k] = 0.;
      }
    }
    }*/

  std::cout << "initialize:: SetPlatePotentials" << std::endl;
  ///////////////////////////////////////////////////
  // Set fixed potentials///
  ///////////////////////////////////////////////////
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      V[i][j][0] = V1-500;
      V[i][j][Plate1Position] = V1;
      V[i][j][Plate2Position] = V2;
      V[i][j][nz-1] = 50;
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
      epsilon[i][j][Plate1Position] = epsilon_copper;
      epsilon[i][j][Plate2Position] = epsilon_copper;
      epsilon[i][j][nz-1] = epsilon_copper;
      epsilon[i][j][0]=epsilon_copper;
    }
  }
    
  std::cout << "initialize:: SetHolesPermittivity" << std::endl;
  
  /////////////////////////////////////////////////// 
  // Define other holes in the Teflon layer (vacuum region)
  /////////////////////////////////////////////////// 
  for(int l=0;l<1;l++){
    for(int m=0;m<1;m++){
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
	  // Skip fixed potential elements

	  if(epsilon[i][j][k] == epsilon_copper) continue;
	  
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



void solveLaplaceEquation()
{

  std::cout << "Solving Laplace" <<  std::endl;
  
  // Create a temporary 3D array for updated potentials
  double*** Vnew = new double**[nx];
  for (int i = 0; i < nx; ++i) {
    Vnew[i] = new double*[ny];
    for (int j = 0; j < ny; ++j) {
      Vnew[i][j] = new double[nz];
      for (int k = 0; k < nz; ++k) {
	Vnew[i][j][k] = V[i][j][k]; // initialize with current potential
      }
    }
  }
  
  int iter = 0;
  double maxDiff = 0.0;

  // Iterative relaxation loop
  do {
    maxDiff = 0.0;
    // Loop over all interior points (assuming boundaries are handled separately)
    for (int k = 1; k < nz - 1; ++k) {
      
      for (int j = 1; j < ny - 1; ++j) {
	for (int i = 1; i < nx - 1; ++i) {
	  // Skip fixed-potential copper points.

	  if (epsilon[i][j][k] == epsilon_copper) continue;
	  
	  // Weighted contributions from the six neighbors.
	  // Note: Adjust weighting if needed for your discretization.
	  double weightSum = 0.0;
	  double weightedV = 0.0;
	  
	  // X-direction neighbors
	  weightedV += epsilon[i+1][j][k] * V[i+1][j][k];
	  weightSum += epsilon[i+1][j][k];
	  weightedV += epsilon[i-1][j][k] * V[i-1][j][k];
	  weightSum += epsilon[i-1][j][k];
	  
	  // Y-direction neighbors
	  weightedV += epsilon[i][j+1][k] * V[i][j+1][k];
	  weightSum += epsilon[i][j+1][k];
	  weightedV += epsilon[i][j-1][k] * V[i][j-1][k];
	  weightSum += epsilon[i][j-1][k];
	  
	  // Z-direction neighbors
	  weightedV += epsilon[i][j][k+1] * V[i][j][k+1];
	  weightSum += epsilon[i][j][k+1];
	  weightedV += epsilon[i][j][k-1] * V[i][j][k-1];
	  weightSum += epsilon[i][j][k-1];
	  
	  double newV = weightedV / weightSum;
	  maxDiff = std::max(maxDiff, std::fabs(newV - V[i][j][k]));
	  Vnew[i][j][k] = newV;
	}
      }
    }
    
    // Copy new potentials back into V for the next iteration.
    for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
	for (int k = 1; k < nz - 1; ++k) {
	  // Only update non-fixed points.
	  if (epsilon[i][j][k] != epsilon_copper) {
	    V[i][j][k] = Vnew[i][j][k];
	  }
	}
      }
    }
    
    iter++;
    // Optionally, print iteration progress
    std::cout << "Iteration " << iter << ", max difference = " << maxDiff << std::endl;
  } while (maxDiff > tolerance && iter < max_iterations);
  
  // Clean up temporary array
  for (int i = 0; i < nx; ++i) {
    for (int j = 0; j < ny; ++j) {
      delete[] Vnew[i][j];
    }
    delete[] Vnew[i];
  }
  delete[] Vnew;
  
  std::cout << "Converged in " << iter << " iterations with max diff " << maxDiff << std::endl;
}





void solveLaplaceSOR()
{
    // Over-relaxation parameter, typically 1 < omega < 2. Adjust for optimal performance.
    double omega = 1.7;
    int iter = 0;
    double max_diff;

    do {
      max_diff = 0.0;
      // Loop over interior points. Adjust if your boundaries need special handling.
      for (int k = 1; k < nz - 1; ++k) {   // Move k-loop to outermost
	for (int j = 1; j < ny - 1; ++j) {
	  for (int i = 1; i < nx - 1; ++i) {
	    // Skip fixed copper points.
	    if (std::abs(epsilon[i][j][k] - epsilon_copper) < 1e-6)  continue;
            
	    // Calculate weighted contributions from the six neighbors.
	    // The weighting uses the epsilon of each neighbor.
	    double weightSum = 
	      epsilon[i+1][j][k] + epsilon[i-1][j][k] +
	      epsilon[i][j+1][k] + epsilon[i][j-1][k] +
	      epsilon[i][j][k+1] + epsilon[i][j][k-1];
	    
	    double neighborSum = 
	      epsilon[i+1][j][k] * V[i+1][j][k] +
	      epsilon[i-1][j][k] * V[i-1][j][k] +
	      epsilon[i][j+1][k] * V[i][j+1][k] +
	      epsilon[i][j-1][k] * V[i][j-1][k] +
	      epsilon[i][j][k+1] * V[i][j][k+1] +
	      epsilon[i][j][k-1] * V[i][j][k-1];
	    
	    // Compute the new potential value from the weighted average.
	    double newV = neighborSum / weightSum;
	    double diff = newV - V[i][j][k];
            
	    // Update using SOR.
	    V[i][j][k] += omega * diff;
            
	    // Track the maximum change over the grid.
	    max_diff = std::max(max_diff, std::fabs(omega * diff));
	  }
	}
      }
      
      iter++;
      // Optionally, output iteration status:
      std::cout << "Iteration " << iter << ", max difference = " << max_diff << std::endl;
    } while (max_diff > tolerance && iter < max_iterations);
    
    std::cout << "SOR converged in " << iter << " iterations with max diff " << max_diff << std::endl;
}


// Conjugate Gradient Solver for Laplace Equation
void solveLaplaceCG() {

  double residual_threshold = tolerance;
  
  // Residual, direction, and search vectors
  std::vector<std::vector<std::vector<double>>> r(nx,std::vector<std::vector<double>>(ny, std::vector<double>(nz, 0.0)));
  std::vector<std::vector<std::vector<double>>> p = r;
  std::vector<std::vector<std::vector<double>>> Ap = r;
  
  // Compute initial residual r = b - Ax (where A is the Laplacian operator)
  double residual_norm = 0.0;
  
  for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int k = 1; k < nz - 1; ++k) {
	if (epsilon[i][j][k] == epsilon_copper) continue;
	
	double laplace_V = (V[i+1][j][k] + V[i-1][j][k] +
			    V[i][j+1][k] + V[i][j-1][k] +
			    V[i][j][k+1] + V[i][j][k-1] - 6.0 * V[i][j][k]);
	
	r[i][j][k] = -laplace_V;
	p[i][j][k] = r[i][j][k];
	residual_norm += r[i][j][k] * r[i][j][k];
      }
    }
  }
  
  double residual_init = residual_norm;

  residual_threshold=1e-6;
  
  
  std::cout << "Initial residual_norm: " << residual_norm << std::endl;
  std::cout << "Initial residual_thr: " << residual_threshold << std::endl;
  std::cout << "Initial residual_init: " << residual_init << std::endl;

  int iter = 0;
  while (iter < max_iterations && residual_norm > residual_threshold * residual_init) {
    double alpha, beta, Ap_dot_p = 0.0, r_dot_r = residual_norm;
    
    // Compute A * p
    for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
	for (int k = 1; k < nz - 1; ++k) {
	  if (epsilon[i][j][k] == epsilon_copper) continue;
	  
	  Ap[i][j][k] = (p[i+1][j][k] + p[i-1][j][k] +
			 p[i][j+1][k] + p[i][j-1][k] +
			 p[i][j][k+1] + p[i][j][k-1] - 6.0 * p[i][j][k]);

	  Ap_dot_p += p[i][j][k] * Ap[i][j][k];
	}
      }
    }
    
    // Compute step size alpha
    alpha = r_dot_r / Ap_dot_p;
    
    // Update solution V and residual r
    residual_norm = 0.0;
    for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
	for (int k = 1; k < nz - 1; ++k) {
	  if (epsilon[i][j][k] == epsilon_copper) continue;
	  
	  V[i][j][k] += alpha * p[i][j][k];
	  r[i][j][k] -= alpha * Ap[i][j][k];
	  
	  residual_norm += r[i][j][k] * r[i][j][k];
	}
      }
    }
    
    // Compute beta for next iteration
    beta = residual_norm / r_dot_r;
    
    // Update direction p
    for (int i = 1; i < nx - 1; ++i) {
      for (int j = 1; j < ny - 1; ++j) {
	for (int k = 1; k < nz - 1; ++k) {
	  if (epsilon[i][j][k] == epsilon_copper) continue;
	  
	  p[i][j][k] = r[i][j][k] + beta * p[i][j][k];
	}
      }
    }

    iter++;
    std::cout << "Iteration " << iter << ", Residual: " << residual_norm << std::endl;
  }
  
  std::cout << "CG Converged in " << iter << " iterations.\n";
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

    double Ex,Ey,Ez;
    
    for (int i = 2; i < nx - 2; i += 1) {
      for (int j = 2; j < ny - 2; j += 1) {
	for (int k = 2; k < nz - 2; k += 1) {
	  
	  // Compute electric field components using fourth-order differences
	  Ex = -(-V[i+2][j][k] + 8.0 * V[i+1][j][k] - 8.0 * V[i-1][j][k] + V[i-2][j][k]) / (12.0 * grid_spacing);
	  Ey = -(-V[i][j+2][k] + 8.0 * V[i][j+1][k] - 8.0 * V[i][j-1][k] + V[i][j-2][k]) / (12.0 * grid_spacing);
	  Ez = -(-V[i][j][k+2] + 8.0 * V[i][j][k+1] - 8.0 * V[i][j][k-1] + V[i][j][k-2]) / (12.0 * grid_spacing);
	  
	  double magnitude = sqrt(Ex * Ex + Ey * Ey + Ez * Ez);
	  
	  if(j == ny/2){
	    fmagn->SetBinContent(i,k,magnitude);
	    
	    if (magnitude > 1 && i%3==0 && k%3==0) { // Avoid plotting zero fields
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

  double Ex, Ey, Ez;
  
  // Loop over all i, j, k (excluding boundaries)
  for (int i = 1; i < nx - 1; ++i) {
    for (int j = 1; j < ny - 1; ++j) {
      for (int k = 1; k < nz - 1; ++k) {	
	
	Ex = -(-V[i+2][j][k] + 8.0 * V[i+1][j][k] - 8.0 * V[i-1][j][k] + V[i-2][j][k]) / (12.0 * grid_spacing);
	Ey = -(-V[i][j+2][k] + 8.0 * V[i][j+1][k] - 8.0 * V[i][j-1][k] + V[i][j-2][k]) / (12.0 * grid_spacing);
	Ez = -(-V[i][j][k+2] + 8.0 * V[i][j][k+1] - 8.0 * V[i][j][k-1] + V[i][j][k-2]) / (12.0 * grid_spacing);
	
	// Calculate the position
	double x = (i - nx / 2) * grid_spacing / 10; //in cm
	double y = (j - ny / 2) * grid_spacing / 10; //in cm
	double z = (k - nz / 2) * grid_spacing / 10; //in cm
	
	// Write to the file
	outfile << x << " " << y << " " << z << " " << Ex*10 << " " << Ey*10 << " " << Ez*10 << "\n";
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
  //solveLaplace();
  //solveLaplaceEquation();
  //solveLaplaceSOR();
  solveLaplaceCG();
  //visualizePotential("Final_Potential","Plot_Final_Potential");
  //visualizePotential();
  visualizeElectricFieldXZProjection();
  writeElectricFieldToFile("OutFileField/field_test.txt");
  theApp.Run();
}

