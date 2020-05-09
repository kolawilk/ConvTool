# ConvTool
Conversion of variables of a flow configuration into dimensionless lattice units for use with PALABOS

Developed at University of Applied Sciences Duesseldorf @ Software Lab, Winter Term 2019 / 2020

## Authors
Lars Wilke, lars.wilke@study.hs-duesseldorf.de
Daniel Morton, daniel.morton@study.hs-duesseldorf.de
Tim Kramer, tim.kramer@hs-duesseldorf.de

## Installation instructions
* Clone this repository with submodules:
```
git clone --recursive https://github.com/kolawilk/ConvTool
```

### Build CoolProp
See http://www.coolprop.org/coolprop/wrappers/StaticLibrary/index.html#static-library for further instructions.

Build on Mac e.g.:

First open './CoolProp/CmakeLists.txt' in your texteditor and uncomment the 'DARWIN_USE_LIBCPP' option, save the file. Then proceed:

'''
cd ConvTool/CoolProp/
mkdir -p build && cd build
cmake .. -DCOOLPROP_STATIC_LIBRARY=ON
cmake --build .
'''

CoolProp can be used now.

### Compile ConvTool

In your terminal cd back to the repository root folder if not already done:

'''
cd ../../
'''

The ConvTool can now be compiled using the following compiler flags:

'''
g++ -L ./CoolProp/build/ -stdlib=libc++ -o ConvTool.o ConvTool.cpp -I ./CoolProp/include/ -I ./CoolProp/externals/fmtlib/ -I ./CoolProp/build/ -lCoolProp
'''

Note: Normally 3 warnings were generated, these can be ignored because they have no influence on the conversion.

The ConvTool can now be executed via './ConvTool.o'.
## Instruction for usage
