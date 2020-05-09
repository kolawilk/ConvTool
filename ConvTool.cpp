//
//  Conversion_Tool
//  Reads in config files with flow setup and converts them into dimensionless lattice variables
//  Writes "palabos_" setup file, which can be read into a PALABOS simulation, if readinpalabos.cpp is included in simulation setup
//  See simrum_2d_convtool.cpp for a full example setup
//

#include "CoolProp.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include "sys/stat.h"
#include <errno.h>
#include <dirent.h>

using namespace std;
using namespace CoolProp;

/* **********************
 * *  Global variables  *
 * **********************/
string fluid, project;
float T_p, T_l, p_p, p_l, Re_p, Re_l, u_p, h_p, l_p, rho_p, visc_p, visc_l, c_s_p, Ma_p, delta_x, tau_l, C_l, C_u, delta_t, secondsteps, Ma_l;
int h_l, l_l;
float c_s_l = sqrt(1.0/3.0);
float u_l;
float Re_g;
int cid, out_fp;
fstream config, output;
int readconf(string), fluidcalc(), latticeconv(), writepalabos(string);
vector<string> files = vector<string>();
//end

/*
 Function to convert the physical units into dimensionless lattice units
 Based on stability or accuracy config
 */
int latticeconv() {
    //Init values to get into loop
    u_l = 1.0;
    Re_g = 11.0;
    delta_x = min(h_p,l_p); //yield at least N=1
    float decrement = delta_x/1000;
    // cid == 1 for accuracy
    if (cid == 1) {
        while (u_l > 0.1) {
            if (delta_x-decrement < 0) {
                decrement /= 1000;
            }
            delta_x -= decrement;
            visc_l = pow(c_s_l,2.0)*(tau_l-0.5);
            delta_t = pow(c_s_l,2.0)*(tau_l-0.5)*(pow(delta_x,2.0))/(visc_p);
            secondsteps = 1.0/delta_t;

            C_u = delta_x/delta_t;

            u_l = u_p/C_u;
            Ma_l = u_l/c_s_l;
            Re_l = u_l*h_l/(visc_l);
        }
        C_l = delta_x;
        l_l = l_p/C_l+0.5;
        h_l = h_p/C_l+0.5;
    }
    // cid == 2 for stability
    else if (cid == 2) {
        while (u_l > 0.4 || Re_g > 10) {
            if (delta_x-decrement < 0) {
                decrement /= 1000;
            }
            delta_x -= decrement;
            visc_l = pow(c_s_l,2.0)*(tau_l-0.5);
            delta_t = pow(c_s_l,2.0)*(tau_l-0.5)*(pow(delta_x,2.0))/(visc_p);
            secondsteps = 1.0/delta_t;

            C_u = delta_x/delta_t;

            u_l = u_p/C_u;
            Ma_l = u_l/c_s_l;
            Re_l = u_l*h_l/(visc_l);
            Re_g = u_l*delta_x/(visc_l);
        }
        C_l = delta_x;
        l_l = l_p/C_l;
        h_l = h_p/C_l;
    }
    // ensure at least N=1
    if (h_l<1) {
        h_l = 1;
    }
    if (l_l<1) {
        l_l = 1;
    }
    return 1;
} //eofunc

/*
 Function to write the converted values to a file that can be read in with with PALABOS
 Exports in ./converted/"projectname".txt
 */
int writepalabos(string outfile) {
    output.open(outfile,ios::out);
    if(!output)
    {
      cout<<"Error in creating file.."<<endl;
      return 0;
    }
    output << "#Input File for PALABOS created with Conversion Tool" << endl;
    output << "project=" << project << endl;
    output << "omega=" << 1/tau_l << endl;
    output << "xlattice=" << l_l << endl;
    output << "ylattice=" << h_l << endl;
    output << "umax=" << u_l << endl;
    output << "deltax=" << delta_x << endl;
    output << "deltat=" << delta_t << endl;
    output << "reynolds=" << Re_p << endl;
    output << "gridreynolds=" << Re_g << endl;
    output << "machlat=" << Ma_l << endl;
    output << "# Conversion parameters:" << endl;
    output << "velocityconv=" << C_u << endl;
    output << "xyconv=" << C_l << endl;
    output.close();
    cout << "PALABOS file created succesfully: " << outfile << endl;
    return 1;
    } //eofunc

/*
 Function to read in all variables defined in config.txt file.
 For an example on defined variables see exampleconfig.txt
*/
int readconf(string conffile) {
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
        if (name.find("reynolds") != string::npos) {
            Re_p = ::atof(value.c_str());
        }
        else if (name.find("pressure") != string::npos) {
            p_p = ::atof(value.c_str());
        }
        else if (name.find("temperature") != string::npos) {
            T_p = ::atof(value.c_str());
        }
        else if (name.find("velocity") != string::npos) {
            u_p = ::atof(value.c_str());
        }
        else if (name.find("height") != string::npos) {
            h_p = ::atof(value.c_str());
        }
        else if (name.find("length") != string::npos) {
            l_p = ::atof(value.c_str());
        }
        else if (name.find("fluid") != string::npos) {
            fluid = value;
        }
        else if (name.find("criteria") != string::npos) {
            cid = ::atof(value.c_str());
        }
        else if (name.find("output") != string::npos) {
            out_fp = ::atof(value.c_str());
        }
        else if (name.find("tau") != string::npos) {
            tau_l = ::atof(value.c_str());
        }
        else if (name.find("project") != string::npos) {
            project = value;
        }
    }
    config.close();
    return 1;
} //eofunc

/*
 Function to calculate needed fluid properties with CoolProp library based on input values of config file
 */
int fluidcalc() {
    rho_p = PropsSI("D","T",T_p,"P",p_p,fluid); // density
    visc_p = PropsSI("V","T",T_p,"P",p_p,fluid)/rho_p; //kinematic viscosity
    c_s_p = PropsSI("A","T",T_p,"P",p_p,fluid); //speed of sound
    if (u_p > 0 && Re_p > 0) {
        cout << "Velocity and Reynolds number defined in config. Abort." << endl;
        return 0;
    }
    else if (u_p > 0) {
        Re_p = u_p*h_p/(visc_p);
    }
    else if (Re_p > 0) {
         u_p = Re_p*visc_p/h_p;
    }
    else {
        cout << "No Reynolds number or velocity defined in config. Abort." << endl;
        return 0;
    }
    Ma_p = u_p/c_s_p;
    //Prints calculated fluid variables to console if output defined in config
    if (out_fp==1) {
        cout << "Calculated fluid properties for project " << project << " with CoolProp:" << endl;
        cout << "Velocity: " << u_p << " m/s" << endl;
        cout << "Density: " << rho_p << " kg/m^3" << endl;
        cout << "Kinematic Viscosity: " << visc_p << " m^2/s" << endl;
        cout << "Speed of Sound: " << c_s_p << " m/s" << endl;
        cout << "Mach number: " << Ma_p << " --" << endl;
        cout << "Reynolds: " << Re_p << " --" << endl;
    }
    return 1;
} //eofunc

/*
 Function to get all files in given path
 based on https://stackoverflow.com/questions/24447540/c-how-to-search-files-in-a-directory-with-certain-name
 */
int getdir (string dir, vector<string> &files)
{
    DIR *dp;
    struct dirent *dirp;

    if((dp = opendir(dir.c_str())) == NULL)
    {
      cout << "Error(" << errno << ") opening " << dir << endl;
      return errno;
    }
    while ((dirp = readdir(dp)) != NULL) {
    string fname = dirp->d_name;
    if(fname.find(".txt") != std::string::npos)
        files.push_back(fname);
    }
    closedir(dp);
    return 0;
} //eofunc

int main()
{
    //Declaration of variables
    string conffind = "config";
    string outname;
    string outpath;
    string confpath;
    int abort;
    cout << "Conversion-Tool from physical units to dimensionless lattice units" << endl << "for 2D simulation with Lattice Boltzmann method with a BGK operator." << endl << endl;
    //Get path with config files
    cout << "Paste *PATH* to config.txt file: ";
    cin >> confpath;
    //get all files in given path
    getdir(confpath,files);
    //defines output path to ./converted/
    outpath = confpath+"/converted/";
    //make ./converted/ directory
    mkdir(outpath.c_str(), S_IRUSR | S_IWUSR | S_IXUSR);
    //convert setup of each config file
    for (unsigned int i = 0;i < files.size();i++) {
        //read in config file
        abort = readconf(confpath+"/"+files[i]);
        //calculate needed fluid properties
        abort = fluidcalc();
        //check if readconf oder fluidcalc fails
        if (abort == 0) {
            return 0;
        }
        //convert into lattice units
        latticeconv();
        //define output name
        outname = "/palabos_"+project+".txt";
        //write setup to read in with palabos
        writepalabos(outpath+outname);
    }
    cout << "Conversion done. " << endl;
    return 1;
}

//eof
