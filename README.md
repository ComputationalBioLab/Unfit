# Unfit
Data fitting and optimization software

### Getting Started (User Version)

The *User* version will build Unfit and run a simple example with the different optimisers. You can use this to build Unfit as a library too.

1. Download & install Code::Blocks (http://www.codeblocks.org/).
   Note: if you are on Windows, choose the version with GCC (mingw)
2. Download Unfit, e.g.
   
   ```
   git clone https://github.com/ComputationalBioLab/Unfit.git
   ```

3. Start Code::Blocks and open the Unfit project file (Unfit.cbp)
4. Build & Run

### Getting Started (Testing Version)

The *Testing* version will build Unfit and run all of the unit tests we use to check Unfit along with a decent number of examples that include fitting different types of problems. You can also use this to build Unfit as a library (does not include the tests).

1. Download & install Code::Blocks (http://www.codeblocks.org/).
   Note: if you are on Windows, choose the version with GCC (mingw)
2. Download Unfit & UnitTest++, e.g.

   ```
   git clone https://github.com/ComputationalBioLab/UnitTest-cpp.git
   git clone https://github.com/ComputationalBioLab/Unfit.git
   ```

3. Start Code::Blocks and open the UnitTest++ project file (UnitTest++.cbp)
4. Build (this should create Debug and Release libraries by default).
   Note: you only need to do this the first time you want to compile Unfit
5. Close the UnitTest++ project
6. Open the Unfit test project file (TestUnfit.cbp)
7. Build & Run

### What Now?

Create your own cost functions so you can minimise what is important to you. Even if you don't want to compile the test version, it may be useful to have a look at the examples to get ideas as to how to code up your custom cost functions. 

We have created a [tutorial](https://computationalbiolab.github.io/assets/pdf/unfit_tutorial.pdf) which provides a step by step guide to writing your own cost functions from scratch. In case you are interested, we have also created another [tutorial](https://computationalbiolab.github.io/assets/pdf/unfit_write_an_optimizer_tutorial.pdf) to explain how to write your own optimizer using the Unfit framework. Both of these are available from the Unfit website.

### Development Environment

Unfit has been developed primarily on GNU Linux using various flavours of Ubuntu. For the current release, testing and development was done on Ubuntu. Our IDE of choice is Code::Blocks, and our compiler GCC (and mingw). We also compile with clang and on Windows from time to time to check this work there too. Please note that the core of Unfit requires C++11 support. We use the UnitTest++ testing framework and offer our sincere thanks to the developers of this software...

### UnitTest++

UnitTest++ is distributed with Unfit for your convenience. All of the credit for UnitTest++ should go to the original authors, who allow redistribution of their source. The license for UnitTest++ is different from that of Unfit, and can be found in the COPYING file in the UnitTest++ directory. The code in the UnitTest++ directory was downloaded from sourceforge.net. They have now moved to github and we encourage you to take a look there for any updates:

https://github.com/unittest-cpp/unittest-cpp

Unfit was developed using UnitTest++ version 1.4 without any major modifications. The only thing we did was to create the Code::Blocks project file so that the UnitTest++ library can be easily compiled alongside Unfit using the same IDE. We have, however, written some nice utility functions to allow you to easily run just one test, or one suite. One day we will get around to upgrading. 

### A Note to Clang Users

We do compile Unfit with clang from time to time, so it is likely that the code will compile with clang++. Consider it an unsupported feature. On Linux, we usually use the conversion tool *cbp2make* to convert the CodeBlocks project files into make files. After that, we simply use a text editor to replace gcc with clang/clang++/llvm-ar as needed in the make files, and build using the command line. If we are feeling enthusiastic we also utilise clang/llvm's static analyser to make sure our code is clean.

### A Note to Visual Studio Users

If you are a hard core Visual Studio fan and don't want to use Code::Blocks, it is possible to compile Unfit in Visual Studio. We do it on a very irregular basis. Most of the code works as advertised when compiling with Visual Studio. However, there is one key thing to be aware of and that is the random numbers generated from a uniform integer distribution do not follow the same sequence as they do with Linux/GCC. The Genetic Algorithm method is the worst affected by this. It still works, but some of our tests fail because e.g. a different number of mutations are generated in different chromosomes. The other methods that is affected is Particle Swarm. To reiterate, these methods still work, it is just the test output that is different, and we don't want to add any platform specific compiler directives to the code.

### A Note to Mac Users

As per the above note to Visual Studio users, this random number issue also appears on the Mac. Some tests may fail due to the stochastic nature of some of the algorithms, but the underlying method is still valid. Note that we don't have a Mac platform to test on and so are probably not much use for Mac support.

### Contributing

We are serious about this being around for a long time and hopefully becoming useful to more and more people. After much study, we have decided to require contributing license agreements (CLAs) from our contributors to keep our work open. Our CLA is adapted directly from the one used by the Apache Foundation, and can be found [here](https://computationalbiolab.github.io/assets/cla/IndividualContributorLicenseAgreement.pdf). If you are interested in how a CLA protects you as well as us, a good explanation is given [here](https://julien.ponge.org/blog/in-defense-of-contributor-license-agreements/). In case you can't be bothered reading the whole thing, here is our summary of the key points of the CLA:

1. (Just definitions)
2. You still own what you contribute, and can do **whatever you want** with it, but you allow us to use it forever, so we don't have to worry that you will change your mind and take it back
3. You promise you are not going to contribute something, then patent it and come after us for patent infringement
4. You are allowed to contribute what you are contributing
5. You are submitting something you have written, not something that belongs to someone else
6. You don't have to provide support unless you choose to do so

If you just want to use our work, you don't have to sign anything - just use it. If you want to contribute please sign, scan and email it back to us. You only have to do this once to be able to contribute to any of our projects.


