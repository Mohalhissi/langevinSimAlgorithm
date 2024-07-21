#include<iostream>
#include<fstream>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<iomanip>
#include"gaussianRandomNum.h"
#include"thetaFunction.h"
#include"initializeParticleCoord.h"
#include"computeForce.h"
#include"integrator.h"

using namespace std;

int main(int argc, char *argv[]){

    cout.precision(10);
    cout<<fixed;
    char *ptr;
    long dumm=strtol(argv[1], &ptr, 10);
    // Convert string to long integer. dumm is used to count the different simualtion outputs.

    long idum= -1* dumm ; // idum is necessary for the random number generator


    // Declaring the varibales describing the system phyiscs:

    double x0=0; //particle initial position stored once after equiliration
    double x=0 ; //particle position (dynamical variable)
    double vx0=0; //particle initial velocity stored once after equiliration
    double vx=0; //particle position (dynamical variable)
    double dt=0.001; // timestep size
    double fx=0; // total force acting on the particle due to the potential
    double en=0,  ke=0, etot=0; // potential energy, kinetic energy, total energy
    double randz=0, randth=0; // two Gaussian andom numbers
    double k_B=1; //Boltzmann constant
    double mass=1; // particle mass
    double gamma=20; // medium's drag coefficient
    double Eb=5; // energy barrier height of the potential
    double temperature=1.0 ;
    double K_t=0, autocorrX_t=0, autocorrVx_t=0, MSD=0; // dynamical correlation functions
    double ave=0, avesq=0, summ=0, summsq=0, variance=0;
    int counter10steps = 1;
    int counter100steps = 0;
    int counter1000steps = 0;

    int const maxbin=500 ;
    double deltaHistogram=10.0/maxbin ;
    double histog_x[maxbin]={0} ; // histogram to bin particle positions
    double histog_vx[maxbin]={0}; // histogram to bin particle velocities
    double sigma= sqrt(2*k_B*temperature*gamma/mass); // dynamical variable

    int EquilibrationTime= 6e6 ;
    int totalSimulationTime= 12e6;
    int r=0, m=0;// integers to output data on a logarthmic time scale

    ofstream vxpdf,measure,simulation, potential;

    vxpdf.open("velocityDistribution.dat" , ios::out) ;
    simulation.open("log.dat" , ios::out);
    measure.open("measurementPhase.dat" , ios::out) ;
    potential.open("potentialFunction.dat" , ios::out) ;

    vxpdf.precision(10); vxpdf<<fixed;
    simulation.precision(10); simulation<<fixed;
    measure.precision(10); measure<<fixed;
    potential.precision(10); potential<<fixed;

    const char *Directory="runsOutputData/";
    const char *Filename=".runsOutputData.dat" ;
    char name_buffer[512];
    ofstream runss;
    runss.precision(10); runss<<fixed;
    sprintf(name_buffer, "%s%ld%s", Directory,dumm,Filename);
    runss.open(name_buffer, ios::out);
    runss.precision(10); runss<<fixed;

    init(x,vx,idum) ; // initilize the algorithm by generating a random position and velocity

    // Now the simulation begins when the loop starts
    // (first equilibration phase starts and then measurement phase follows)

    simulation<<"simulation time"<<setw(25)<<"temprature Fluctuation"<<setw(15)<<"position"
        <<setw(15)<<"velocity"<<setw(25)<<"potential energy"<<setw(25)<<"kinetic energy"<<setw(25)<<"total energy"<<endl ;


    for (int timestep=0 ; timestep<=totalSimulationTime ; ++timestep) {

        force(fx,en,x,mass,Eb) ; // compute the force on the particle

        randz=gasdev(idum) ; // two random numbers for the thermal random force computation
        randth=gasdev(idum) ;

        integrate(x,fx,en,vx,randz,randth,gamma,sigma,etot,dt,mass,Eb) ; // update the position and velocity

        summ=summ+vx ;
        summsq=summsq + pow(vx,2) ;
        ave=summ/(timestep) ;
        avesq=summsq/(timestep) ;
        variance = avesq - pow(ave,2) ;
        double TeFluctutation=mass*variance ; // observe the temperature fluctuations

        ke=0.5*pow(vx,2) ;

        //output simualtion data every 10 timesteps
        if ( timestep/10.0 == counter10steps ){

            cout<< timestep*dt<<setw(15)<<TeFluctutation<<setw(15)<<x
                <<setw(15)<<vx<<setw(15)<<en<<setw(15)<<ke<<setw(15)<<etot<<endl ;

            simulation<<timestep*dt<<setw(15)<<TeFluctutation<<setw(15)<<x<<setw(15)
                      <<vx<<setw(15)<<en<<setw(15)<<ke<<setw(15)<<etot<<endl ;

            counter10steps++;
        }

        /* Check if the equilibration time has been enough or not.
           Otherwise the equilibration time should be increased. If the system is
           equilibrated, the measurement phase directly begins as below*/

        if (timestep==EquilibrationTime){
            x0=x ; vx0= vx ; // recodring the initial coordinates in the measurement phase
        }

        // Now the simulation is the measurement simulation

        if (timestep >= EquilibrationTime && timestep<= totalSimulationTime){

            int ibinx=round(x/deltaHistogram);
            if (ibinx >= -250 && ibinx <= 250) {
                int cc=ibinx+ 250 ;
                histog_x[cc]=histog_x[cc] + 1 ;
            }

            int ibinvx=round(vx/deltaHistogram);
            if (ibinvx >= -250 && ibinvx <= 250) {
                int dd=ibinvx+ 250 ;
                histog_vx[dd]= histog_vx[dd] + 1 ;
            }


            MSD=pow((x-x0),2) ;

            K_t = 4*(theta(x0) - 0.5)*(theta(x)-0.5);

            autocorrX_t=x*x0/(x0*x0);
            autocorrVx_t=vx*vx0/(vx0*vx0);


            if ((timestep-EquilibrationTime) == 0 ){r=0; m=0;}


            if ( (timestep-EquilibrationTime)/100.0 == counter100steps ){
                potential<< x <<setw(15)<<en<<setw(15)<<fx<<endl;
                counter100steps++;

            }

            if ( (timestep-EquilibrationTime)/1000.0 == counter1000steps ){

                measure << (timestep-EquilibrationTime)*dt <<setw(15)
                        <<TeFluctutation<<setw(15)<<x<<setw(15)<<MSD<<setw(15)<<vx<<setw(15)<<fx \
                        <<setw(15)<<en<<setw(15)<<ke<<setw(15)<<etot<<endl ;
                counter1000steps++;

            }

            if ( (timestep-EquilibrationTime) <= 9  && (timestep-EquilibrationTime) <= r*pow(10,m) ) {
                runss <<(timestep-EquilibrationTime)*dt <<setw(15)<<MSD<<setw(15)<< K_t\
                     <<setw(15)<<autocorrX_t<<setw(15)<< autocorrVx_t <<setw(15)<< TeFluctutation <<endl;
                r=r+1;
                if (r==10) {
                    r=1;
                    m=m+1;
                }

            }

            if(timestep-EquilibrationTime == 10) {r=5 ; m=0;}


            if ( (timestep-EquilibrationTime) >= 10  && (timestep-EquilibrationTime)==2*r*pow(10,m) ) {


                runss <<(timestep-EquilibrationTime)*dt<<setw(15)<<MSD<<setw(15)<< K_t\
                     <<setw(15)<<autocorrX_t<< setw(15)<< autocorrVx_t <<setw(15)<< TeFluctutation <<endl;

                r=r+1;
                if (r==50) {
                    r=5;
                    m=m+1;

                }

            }

        }

    }

    const char *Directory1="probability/";
    const char *Filename1=".positionProbabilitydistribution.dat" ;
    char name1_buffer[512];

    ofstream xpdf;
    xpdf.precision(10); xpdf<<fixed;
    sprintf(name1_buffer, "%s%ld%s", Directory1,dumm,Filename1);
    xpdf.open(name1_buffer, ios::out);

    for (int ibinx=0 ; ibinx<500 ; ++ibinx){

        int cc=ibinx-250 ;

        xpdf<<cc*deltaHistogram<<setw(15)<<histog_x[ibinx]/(deltaHistogram*(totalSimulationTime-EquilibrationTime))<<endl ;
    }

    for (int ibinvx=0 ; ibinvx<500 ; ++ibinvx){

        int dd=ibinvx-250 ;

        vxpdf<<dd*deltaHistogram<<setw(15)<<histog_vx[ibinvx]/(deltaHistogram*(totalSimulationTime-EquilibrationTime))<<endl ;
    }




    cout << endl;
    cout << "Energy barrier height: Eb= " << Eb << endl;
    cout << "Fluid dag: gamma= " << gamma << endl;
    cout << "timestep size dt="<< dt << endl;
    printf("\nThe simulation time is %d dt \n", totalSimulationTime );
    printf("The equilibration time is %d dt \n", EquilibrationTime );
    printf("The measuring time is %d dt \n", totalSimulationTime-EquilibrationTime);


    simulation.close(); measure.close();
    vxpdf.close(); potential.close(); measure.close();
    runss.close(); xpdf.close();


    cout << endl;
    return(0) ;

}








