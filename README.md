# Unfit
Data fitting and optimization software

### Getting Started (User Version)

1. Download & install Code::Blocks (http://www.codeblocks.org/).
   Note: if you are on Windows, choose the version with GCC (mingw)
2. Download Unfit
    ```git clone https://github.com/ComputationalBioLab/Unfit.git```
3. Start Code::Blocks and open the Unfit project file (Unfit.cbp)
4. Build & Run

### Getting Started (Development Version)

1. Download & install Code::Blocks (http://www.codeblocks.org/).
   Note: if you are on Windows, choose the version with GCC (mingw)
2. Download Unfit
    ```git clone https://github.com/ComputationalBioLab/Unfit.git```
3. Start Code::Blocks and open the UnitTest++ project file (UnitTest++.cbp)
4. Build (this should create Debug and Release libraries by default).
   Note: you only need to do this the first time you want to compile Unfit
5. Close the UnitTest++ project
6. Open the Unfit development project file (Unfit-devel.cbp)
7. Build & Run

### What Now?

Create your own cost functions so you can minimise what is important to you. Even if you don't want to compile the development version, it may be useful to have a look at the test examples to get ideas as to how to code up your custom cost functions. 

We have created a new tutorial which provides a step by step guide to writing your own cost functions from scratch. We have also created another new tutorial to explain how to write your own optimizer using the Unfit framework. Both of these are available from the Unfit website.

### Development Environment

Unfit has been developed primarily on GNU Linux using various flavours of Ubuntu. For the current release, testing and development was done on Ubuntu. Our IDE of choice is Code::Blocks, and our compiler GCC (and mingw). We also compile with clang and on Windows from time to time to check this work there too. Please note that the core of Unfit requires C++11 support. We use the UnitTest++ testing framework and offer our sincere thanks to the developers of this software...

### UnitTest++

UnitTest++ is distributed with Unfit for your convenience. All of the credit for UnitTest++ should go to the original authors, who allow redistribution of their source. The license for UnitTest++ is different from that of Unfit, and can be found in the COPYING file in the UnitTest++ directory. The code in the UnitTest++ directory was downloaded from sourceforge.net. They have now moved to github and we encourage you to take a look there for any updates:

https://github.com/unittest-cpp/unittest-cpp

Unfit was developed using UnitTest++ version 1.4 without any major modifications. The only thing we did was to create the Code::Blocks project file so that the UnitTest++ library can be easily compiled alongside Unfit using the same IDE. We have, however, written some nice utility functions to allow you to easily run just one test, or one suite. One day we will get around to upgrading. 

### A Note to Clang Users

We do compile Unfit with clang from time to time, so it is likely that the code will compile with clang++. Consider it an unsupported feature. On Linux, we usually use the conversion tool *cbp2make* to convert the CodeBlocks project files into make files. After that, we simply use a text editor to replace gcc with clang/clang++/llvm-ar as needed in the make files, and build using the command line. If we are feeling enthusiastic we also utilise clang/llvm's static analyser to make sure our code is clean.

### A Note to Visual Studio Users

If you are a hard core Visual Studio fan and don't want to use Code::Blocks, it is possible to compile Unfit in Visual Studio. However, at the time of writing Visual Studio does not support some C++11 constexpr functionality (which we use). All of the _problem_ code is in Options.hpp and Options.cpp. All you need to do is explicitly put the defaults into the constructor and the Reset methods,
then you can removed the constexpr lines and all should hopefully compile. 

### A Note to non-Linux Users
There appear to be some differences in the random number generation implementations across different platforms. All of the test examples work well on Linux, but we know that some do not on the Mac which appears to be due to differences with the numbers coming out of the uniform distribution. Some examples may fail due to the stochastic nature of some of the algorithms, and we don't have a Mac platform to test on.

