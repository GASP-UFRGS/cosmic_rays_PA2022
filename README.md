# cosmic_rays_PA2022

Geant4 simulation of cosmic showers at the Earth's atmosphere for the UFRGS Portas Abertas 2022. 


## Pre-requisites

* Geant4. Further information about the installation of Geant4 can be found in the [Geant4 Installation Guide (English)](https://mirror.umd.edu/gentoo/distfiles/Geant4InstallationGuide-4.10.7.pdf), and a [tutorial in portuguese](http://lief.if.ufrgs.br/~marcosderos/miscellaneous/Guia_Geant4.pdf).

## Instructions

Inside the project folder, create a `builddir` folder to build the application:

`$ mkdir builddir && cd builddir`

Load the Geant4 environment variables:

`$ source <path>/bin/geant4.sh`

Compile:

`$ cmake ../ && make`

Run the application:

`./cosmic_rays`

Inside the Visualization GUI, you can run a test macro:

`/control/execute comandos_shower`



