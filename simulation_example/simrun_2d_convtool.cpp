#include "palabos2D.h"
#ifndef PLB_PRECOMPILED // Unless precompiled version is used,
  #include "palabos2D.hh"   // include full template code
#endif
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include "sys/stat.h"
#include <errno.h>
#include <dirent.h>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

T omega, umax, deltax, deltat, reynolds, gridreynolds, machlat, velocityconv, xyconv;
plint xlattice, ylattice;
string project;

void domainSetup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                  plint nx, plint ny, T u,
                  OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition )
{
    // Definition Inlet, Outlet und Wände
    Box2D topWall(0, nx-1, ny-1, ny-1);
    Box2D bottomWall (0, nx-1, 0, 0);
    Box2D inlet (0, 0, 1 ,ny-2);
    Box2D outlet (nx-1, nx-1, 1, ny-2);

    // Definition der Art der Randbedingungen
    boundaryCondition.setVelocityConditionOnBlockBoundaries (lattice, topWall );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (lattice, bottomWall );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (lattice, inlet );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (lattice, outlet, boundary::outflow );

    // Werte für Randbedingungen setzen
    setBoundaryVelocity(lattice, topWall, Array<T,2>((T)0.,(T)0.) );
    setBoundaryVelocity(lattice, bottomWall, Array<T,2>((T)0.,(T)0.) );
    setBoundaryVelocity(lattice, inlet, Array<T,2>(u,(T)0.));

    initializeAtEquilibrium (lattice, lattice.getBoundingBox(),(T)1., Array<T,2>(u,(T)0.)  );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGif(BlockLatticeT& lattice, plint iter)
{
    const plint imSize = 600;

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("uNorm", iter, 6),
                               *computeVelocityNorm(lattice),
                               imSize, imSize );
    imageWriter.writeScaledGif(createFileName("logUnorm", iter, 6),
                               *computeLog(*add((T)1.e-8,*computeVelocityNorm(lattice))),
                               imSize, imSize );
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              T dx, T dt, plint iter)
{

    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}

int readconf(string conffile) {
    fstream config;
    string line, key, comment = "#";
    config.open(conffile,ios::in);
    if(!config) {
        cout << "Error loading config file.." << endl;
        return 0;
        
    }
    cout << "Config file loaded succesfully.." << endl << "Read in configuration.." << endl;
    while (!config.eof()) {
        getline(config,line);
        if (line[0] == '#' || line.empty()) {
            continue;
        }
        auto delimiterPos = line.find("=");
        auto name = line.substr(0, delimiterPos);
        auto value = line.substr(delimiterPos + 1);
        if (name.find("omega") != string::npos) {
            omega = ::atof(value.c_str());
        }
        else if (name.find("xlattice") != string::npos) {
            xlattice = ::atof(value.c_str());
        }
        else if (name.find("ylattice") != string::npos) {
            ylattice = ::atof(value.c_str());
        }
        else if (name.find("umax") != string::npos) {
            umax = ::atof(value.c_str());
        }
        else if (name.find("deltax") != string::npos) {
            deltax = ::atof(value.c_str());
        }
        else if (name.find("deltat") != string::npos) {
            deltat = ::atof(value.c_str());
        }
        else if (name.find("reynolds") != string::npos) {
            reynolds = ::atof(value.c_str());
        }
        else if (name.find("gridreynolds") != string::npos) {
            gridreynolds = ::atof(value.c_str());
        }
        else if (name.find("machlat") != string::npos) {
            machlat = ::atof(value.c_str());
        }
        else if (name.find("velocityconv") != string::npos) {
            velocityconv = ::atof(value.c_str());
        }
        else if (name.find("xyconv") != string::npos) {
            xyconv = ::atof(value.c_str());
        }
        else if (name.find("project") != string::npos) {
            project = value;
        }
    }
    config.close();
    return 1;
}

int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);
    string conffile, outfile;
    fstream output;
    output.open(outfile,ios::out);
    pcout << "Paste path to converted *FILE*:" << endl;
    cin >> conffile;
    readconf(conffile);
    string outpath = "./"+project+"/";
    outfile = outpath+"log.txt";
    mkdir(outpath.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    global::directories().setOutputDir(outpath);
    output.open(outfile,ios::out);
    if(!output) {
        cout << "Error in creating file.." << endl;
        return 0;
    }
    pcout<<"<<Running simulation with configuration: " << conffile << ">>\n" << endl;
    output <<"Running simulation with configuration: " << conffile << ">>\n" << endl;
    
    //More settings:
    
    T imSave = (T)1.0;
    if (deltat > imSave) {
        imSave = 2*deltat;
    }
    const T vtkSave = imSave;
    const T maxT = (T)(xlattice/(umax/2)+1)*deltat;
    const T logT = imSave;
    
    //T u = 0.11111111;

    //const T tau         = 0.55;
    //const T deltaT      = 0.00016667;
    //const T deltaX      = 0.0001;
    //const plint Nx      = 2000;
    //const plint Ny      = 100;

    //const T logT     = (T)0.05;
    //const T imSave   = (T)10.0;
    //const T vtkSave  = (T)0.05;
    //const T maxT     = (T)10.0;

    //writeLogFile(parameters, "2D cavity");

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
              xlattice, ylattice,
              new BGKdynamics<T,DESCRIPTOR>(omega) );

    // Welches wählen? Unterschied??
    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
        //boundaryCondition = createInterpBoundaryCondition2D<T,DESCRIPTOR>();


    domainSetup(lattice, xlattice, ylattice, umax/2, *boundaryCondition);

    T previousIterationTime = T();
    T completeSimTime = T();
    global::timer("simTime").restart();
    writeGif(lattice, 0);
    writeVTK(lattice, deltax, deltat, 0);
    
    // Main loop over time iterations.
    for (plint iT=0; iT*deltat<maxT; ++iT) {
        plint nStep;
        global::timer("mainLoop").restart();
        /* uncomment for gif output:
        nStep = imSave/deltat + 0.5;
        if (iT%(nStep)==0 && iT>0) {
            pcout << "Saving Gif ..." << endl;
            writeGif(lattice, iT);
            pcout << endl;
        }
        */
        nStep = vtkSave/deltat + 0.5;
        if (iT%(nStep)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            output << "Saving VTK file ..." << endl;
            writeVTK(lattice, deltax, deltat, iT);
        }

        //nStep = logT/deltat + 0.5;
        //if (iT%(nStep)==0) {
            pcout << "step " << iT
            << "; t=" << iT*deltat;
            output << "step " << iT
            << "; t=" << iT*deltat;
        //}

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
        //nStep = logT/deltat + 0.5;
        //if (iT%(nStep)==0) {
            pcout << "; av energy="
                  << setprecision(10) << getStoredAverageEnergy(lattice)
                  << "; av rho="
                  << getStoredAverageDensity(lattice) << endl;
            output << "; av energy="
            << setprecision(10) << getStoredAverageEnergy(lattice)
            << "; av rho="
            << getStoredAverageDensity(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
            output << "Time spent during previous iteration: "
            << previousIterationTime << endl;
        //}

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
    completeSimTime = global::timer("simTime").stop();
    pcout << "Time for Simulation: " << completeSimTime << endl;
    output << "Time for Simulation: " << completeSimTime << endl;
}
