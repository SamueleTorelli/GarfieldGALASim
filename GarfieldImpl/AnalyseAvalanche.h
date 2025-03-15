#ifndef ANALYSEAVALANCHE_H
#define ANALYSEAVALANCHE_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <Garfield/AvalancheMicroscopic.hh> 

using namespace Garfield;

void AnalyseAval(AvalancheMicroscopic aval,ComponentGrid* field){
  
  const int np = aval.GetNumberOfElectronEndpoints();
  
  //std::cout << "Number of electrons produced in avalanche: " << np << std::endl;
  double x1, y1, z1, t1, e1;
  double x2, y2, z2, t2, e2;
  // Loop over the electrons produced [should be only one!]
  
  for (int ie = 0; ie < np; ie++) {
    
      int status;
      aval.GetElectronEndpoint(ie, x1, y1, z1, t1, e1, x2, y2, z2, t2, e2, status);
      
      std::cout << "Electron" << ie <<" ends at (x, y, z) = ("
		<< x2 << ", " << y2 << ", " << z2
		<< ") with an energy of " << e2 << " eV.\n"
		<< "Status: " << status << "\n";
      std::cout << "Field at the point: " << field->ElectricField(x2, y2, z2)[0] << " "
		<< field->ElectricField(x2, y2, z2)[1] << " "
		<< field->ElectricField(x2, y2, z2)[2] << "\n"; 
      
      if (status == -5) {	
	std::cout << "Electron out of drif medium" << std::endl;
      }
      else if(status == -1){
	std::cout << "Electron out of the sensor volume" << std::endl;
      }
      else {
	std::cout << "Electron " << ie << " ended with a strange status (" << status << ")\n";
	//        << "(x1, y1, z1) = (" << x1 << ", " << y1 << ", " << z1
	//        << "), t1 = " << t1 << ", e1 = " << e1 << "\n"
	//        << "(x2, y2, z2) = (" << x2 << ", " << y2 << ", " << z2
	//        << "), t2 = " << t2 << ", e2 = " << e2 << "\n";
      }  
  }
  

}

#endif // ANALYSEAVALANCHE_H   
