# ConvTool
Conversion of variables of a flow configuration into dimensionless lattice units for use with PALABOS

Developed at University of Applied Sciences Duesseldorf @ Software Lab, Winter Term 2019 / 2020

## Authors
**Lars Wilke** lars.wilke@study.hs-duesseldorf.de

**Daniel Morton** daniel.morton@study.hs-duesseldorf.de

**Tim Kramer** tim.kramer@hs-duesseldorf.de

# Installation instructions
* Clone this repository with submodules:
```
git clone --recursive https://github.com/kolawilk/ConvTool
```

## Build CoolProp
See http://www.coolprop.org/coolprop/wrappers/StaticLibrary/index.html#static-library for further instructions.

### Build on Mac e.g.:

* First open `./CoolProp/CmakeLists.txt` in your texteditor and uncomment the `DARWIN_USE_LIBCPP` option, save the file. Then proceed:

```
cd ConvTool/CoolProp/
mkdir -p build && cd build
cmake .. -DCOOLPROP_STATIC_LIBRARY=ON
cmake --build .
```

* CoolProp can be used now

## Compile ConvTool

* In your terminal cd back to the repository root folder if not already done:

```
cd ../../
```

* The ConvTool can now be compiled using the following compiler flags:

```
g++ -L ./CoolProp/build/ -stdlib=libc++ -o ConvTool.o ConvTool.cpp -I ./CoolProp/include/ -I ./CoolProp/externals/fmtlib/ -I ./CoolProp/build/ -lCoolProp
```

**Note:** *Normally 3 warnings were generated, these can be ignored because they have no influence on the conversion.*

* The ConvTool can now be executed via `./ConvTool.o`.

## Instructions for usage

This conversion tool converts a setup of flow parameters into dimensionless lattice units. Based on a config file the conversion tool reads in the flow parameters, calculates fluid properties with the help of CoolProp and after conversion provides a `palabos_project.txt` that can be read in to simulate the flow with palabos.

The conversion is mainly based on the statements in "The Lattice Boltzmann Method" [1]. In the current version the conversion tool can be used only in the scope of 2D flows using the BGK operator.

# Setup a config file

The flow parameters are set up in different config files. The tool is capable to convert multiple configurations located in one path. The path can be located anywhere on your file system. An example config file is part of this repository:

```
# Config-File for LBM Conversion Tool
# Variable names are not permitted for changes
#
# Projektname
project=WaterStd_Re1000_Acc
# Fluid string:
# For possible fluid strings see CoolProp Documentation:
# http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids
fluid=Water
# Physical Temperature (K):
temperature=293.15
# Physical pressure (Pa):
pressure=101325
# ONE of the following two options (Re or velocity) must be given
# Reynolds Number (=-1 if not known):
reynolds=1000
# Velocity of fluid (m/s) (=-1 if not known):
velocity=-1
# Channel height (m):
height=0.1
# Channel length (m):
length=2.0
# Criteria for conversion
# (1) Accuracy (Umax < 0.1)
# (2) Stability (Umax < 0.4)
criteria=1;
# Choose tau
# Recommended for start:
# Accuracy: tau = 0.9
# Stability: tau = 0.55
tau=0.9
# Output of calculated fluid properties
output=1;
```

**Note:** *The output file name is based on the definition of the project string variable. So make sure to define a different project name for each config file, otherwise the output files will be overwritten by the conversion tool.*

# Conversion with the ConvTool

To perform the conversion you simply need to execute the ConvTool and paste the path where the config files are located. The ConvTool will generate a new directory in this config path `./converted/` with the output files for use in PALABOS.

**Note:** *The ConvTool can help you to set up a whole set of different flow configurations. Nevertheless, a basic understanding of fluid dynamics and LBM is required. To get an understanding of choosing input parameters e.g.* `tau` *and check the output variables we recommend to read chapter 7 in [1].*

# Simulation with PALABOS

PALABOS is included as submodule in this repository. For further information please read the PALABOS project website [2]. A sample simulation is also part of this repository. For your simulations the readconf() function needs to be added to your simulation code:

* Extract from example simulation with readconf() function:

```
[...]
T omega, umax, deltax, deltat, reynolds, gridreynolds, machlat, velocityconv, xyconv;
plint xlattice, ylattice;
string project;
[...]
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
```

The variables of this input are further implemented in the main section of the simulation file, e.g. project name as output folder for "vti"-files. Please see the example simulation file for further information.

## References

[1] Krüger, T.; Kusumaatmaja, H.; Kuzmin, A.; Shardt, O.; Silva, G. & Viggen, E. M.; The Lattice Boltzmann Method; Springer International Publishing, 2017

[2] PALABOS project website: https://palabos.unige.ch


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
