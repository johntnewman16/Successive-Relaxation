/*
AEP 4380 Final Project

Semiconductor Device Modeling

Jack Newman December 8, 2015

Last Edited: Decemeber 13, 2015

Intel(R) Core(TM) i7-4810MQ CPU @ 2.80 GHz

MinGW C++ compiler, GCC 4.8.1

OS: Windows 8.1

Note: Comments are either inline or below the lines they describe
*/
#include <cstdlib> // importing plain c
#include <cmath> // math library

#include <iostream> // stream IO
#include <fstream> // stream file IO
#include <iomanip> // to format the output
// define the following symbol to enable bounds checking
//
#define ARRAYT_BOUNDS_CHECK
#include "arrayt.hpp" //used with permission by Professor Earl Kirkland

int Nx = 320, Ny = 320; //End coordinates for boundaries

double D1 = 10E18, D2 = 10E18, D3 = -10E16;//Doping concentrations in cm^-3
double V_T = 0.026; //Thermal Voltage
double del = 5E-8;//Dioxide Thickness in m
double BigQ = 0; //Trapped Charges at the interface
double es = 11.68, ed = 3.9; //relative permittivity of silicon and silicon dioxide
//double q = 4.8032E-10; //charge of an electron in esu
double q = 1.60217662E-19;//charge of an electron in Couloumbs
double ni = 1.45E10; // intrinsic carrier concentration of silicon
double B = q*ni/es; 
double steps = 0; //counter for use in relax method
double Vg = 0.2; //Applied Voltage at Gate

double initialgreater(double D){
	//helper function to initialize values on grid that have a doping
	//concentration greater than ni
	double phi = V_T*log(D/(2*ni) + sqrt(1 + D*D/(4*ni*ni)));
	return phi;
}

double initialless(double D){
	//helper function to initialize values on grid that have a doping
	//concentration less than ni
	double phi = V_T*log(1/(sqrt(1 + D*D/(4*ni*ni)) - D/(2*ni)));
	return phi;
}


void relax(arrayt<char>& flag, arrayt<double>& grid){
	//Relaxation method, which takes in a flag and grid array of equal
	//of equal dimensions. Will run until relaxation tolerance is reached.
	double w = 1.86; //relaxation parameter
	double diffmax = 0;//maximum change in potential for a given iteration
	double max = 0;//maximum potential on the grid
	double fd;//finite difference value
	double relaxtolerance;//relaxation tolerance
	int i, j, Nmax = 100000;
	//counters, and maximum steps to catch calls which fail to converge


	for(i = 0; i < (Nx + 1); i++){
		for(j = 0; j < (Ny + 1); j++){
			if(flag(i,j) == '0'){
				//finite difference equation on dioxide interface
				fd = (grid(i, j-1)*del*es/ed+del*BigQ/ed+Vg)/(1+del*es/ed);
			}

			else if(flag(i,j) == '1'){
				//finite difference equation on ohmic contact 1
				fd = V_T*log((D1/2.0 + sqrt(D1*D1/4.0 + ni*ni))/ni);
			}

			else if(flag(i,j) == '2'){
				//finite difference equation on ohmic contact 2
				fd = V_T*log((D2/2.0 + sqrt(D2*D2/4.0 + ni*ni))/ni);
			}

			else if(flag(i,j) == '3'){
				//finite difference equation on ohmic contact 3
				fd = V_T*log((D3/2.0 + sqrt(D3*D3/4.0 + ni*ni))/ni);
			}

			else if(flag(i,j) == '4'){
				//finite difference equation on left artificial zone boundary
				fd = grid(i+1, j);
			}

			else if(flag(i,j) == '5'){
				//finite difference equation on right artificial zone boundary
				fd = grid(i-1, j);
			}

			else if(flag(i,j) == 's'){
				//Application of Newton's Method for finite difference
				//equation in the interior source region
				double F = grid(i+1, j) + grid(i-1, j) + grid(i, j+1) + 
					grid(i, j-1) + q*D1/es;
				double psi = 4*grid(i,j)+2*B*sinh(grid(i,j)/V_T) - F;
				double delpsi = 4+2*B/V_T*cosh(grid(i,j)/V_T);
				double chi = -psi/delpsi;
				double phinew = grid(i, j) + chi;
				fd = phinew;
			}
			else if(flag(i,j) == 'd'){
				//Application of Newton's Method for finite difference
				//equation in the interior drain region
				double F = grid(i+1, j) + grid(i-1, j) + grid(i, j+1) + 
					grid(i, j-1) + q*D2/es;
				double psi = 4*grid(i,j)+2*B*sinh(grid(i,j)/V_T) - F;
				double delpsi = 4+2*B/V_T*cosh(grid(i,j)/V_T);
				double chi = -psi/delpsi;
				double phinew = grid(i, j) + chi;
				fd = phinew;
			}

			else if(flag(i,j) == '7'){
				//Application of Newton's Method for finite difference
				//equation in the interior hole carrier region
				double F = grid(i+1, j) + grid(i-1, j) + grid(i, j+1) + 
					grid(i, j-1) + q*D3/es;
				double psi = 4*grid(i,j)+2*B*sinh(grid(i,j)/V_T) - F;
				double delpsi = 4+2*B/V_T*cosh(grid(i,j)/V_T);
				double chi = -psi/delpsi;

				double phinew = grid(i, j) + chi;
				fd = phinew;
			}

			if(abs(grid(i,j) - fd) > diffmax){
				//Obtains the maximum change for a given iteration
				diffmax = abs(grid(i,j) - fd);
			}

			if(abs(grid(i,j)) > max){
				//Obtains the maximum potential in a given iteration
				max = abs(grid(i,j));
			}
			//Relaxation Equation
			grid(i,j) = grid(i,j) + w*(fd - grid(i,j));
		}
	}
	//cout << "max = " << max << endl;
	if(steps > Nmax){
		//if the number of steps becomes greater than NMAX, throw error
		cout << "Max Step Error" << endl;
		exit(0);
	}
	relaxtolerance = 1/1000.0; //relaxation tolerance
	//cout << "relaxtolerance = " << relaxtolerance << setw(15) << "diffmax = "
	//	<< diffmax << endl;
	if(diffmax > relaxtolerance){
		//If the maximum potential change is greater than the relaxation
		//tolerance parameter, add one to steps and run relax again
		steps += 1;
		relax(flag, grid);
	}
	if(diffmax < relaxtolerance){
		//End relax call once diffmax is less than the 
		// relaxation tolerance parameter
		cout << "DONE" << endl;
		cout << w << setw(15) << steps << endl;
		//print relaxation parameter and steps
	}
}

int main(){
	
	ofstream fp;
	ofstream gp;
	
	fp.open( "applycontour_new.dat" ); //open new file for contour data
	
	if( fp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		// catches fp.open failure, and returns appropriate error
		return( EXIT_SUCCESS );
	}

	gp.open( "1dapply_new.dat" ); //open new file for 1d cross-section data
	
	if( gp.fail() ) { // or fp.bad()
		cout << "cannot open file" << endl;
		// catches fp.open failure, and returns appropriate error
		return( EXIT_SUCCESS );
	}


	int i, j; //counters
	arrayt<double> grid(Nx + 1, Ny + 1);//potential mesh
	arrayt<char> flags(Nx + 1, Ny + 1);//flag
	for(j = 0; j < (Ny + 1); j++){
		for(i = 0; i < (Nx + 1); i++){
			if(j == Ny and i >= 60 and i <= (Nx - 60)){
				//flags for dioxide interface
				flags(i, j) = '0';
			}
			else if(j == Ny and i < 60){
				//flag for ohmic contact 1
				flags(i, j) = '1';
			}
			else if(j == Ny and i > (Nx - 60)){
				//flag for ohmic contact 2
				flags(i, j) = '2';
			}
			else if(j == 0){
				//flag for ohmic contact 3
				flags(i, j) = '3';
			}
			else if(i == 0){
				//flag for left artifical boundary
				flags(i, j) = '4';
			}
			else if(i == Nx){
				//flag for right artifical boundary
				flags(i, j) = '5';
			}
			else if(i <= 60 and j >= (Ny - 60)){
				//flag for source region
				flags(i, j) = 's';
			}
			else if(i >= (Nx - 60) and j >= (Ny - 60)){
				//flag for drain region
				flags(i, j) = 'd';
			}
			else{
				//flag for p region
				flags(i, j) = '7';
			}
		}

	}

	for(j = 0; j < (Ny + 1); j++){
		for(i = 0; i < (Nx + 1); i++){
			if(i <= 60 and j >= (Ny-60)){
				//Initialize source region potentials
				grid(i,j) = initialgreater(D1);

			}
			else if(i >= (Nx-60) and j >= (Ny-60)){
				//initialize drain region potentials
				grid(i,j) = initialgreater(D2);

			}
			else{
				//initialize p-region potentials
				grid(i, j) = initialless(D3);

			}
		}
	}
	//Start Relaxation Method
	relax(flags, grid);

	//Save data for contour graph in applycontour
	for(j = Ny; j > 0; j--){//altered for formatting in imshow function
		for(i = 0; i < (Nx + 1); i++){
			fp << grid(i, j) << setw(15);
		}
		fp << endl;
	}
	fp.close();
	
	//Save data for cross-section graph
	for(j = 0; j <= Ny/2; j++){
		gp << j*0.005 << setw(15) //x axis
		<< grid(j, 220) << setw(15) <<//y=1.1
			grid(j, 230) << setw(15) << //y = 1.15
			grid(j, 240) <<//y = 1.2
			setw(15) << grid(j, 250) <<//y=1.25
			setw(15) << grid(j, 258) <<//y=1.29
			setw(15) << grid(j, 265) <<//y=1.325
			setw(15) << grid(j, 295) <<//y=1.475
			setw(15) << grid(j, 315) << endl;//y=1.575
	}
	gp.close();


	return(EXIT_SUCCESS);
}//end main
