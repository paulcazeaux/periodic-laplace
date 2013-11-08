periodic-laplace
================

periodic-laplace is a module for the solution of boundary-element problems in a periodic setting. It uses an approximation to the Green's kernel to the Laplace equation in a cube of size 1 with periodic boundary conditions to write a few boundary operators. The framework used is the open-source boundary element library BEM++, interfaced with M. Bebendorf's AHMED library.

History
=======

0.1 (8 November 2013)
A first developement version using the ideas developed at Cemracs '13 while working on the COMPRESS project, including

- Boundary integral operators (single-layer potential, double-layer potential, adjoint double-layer potential) for the periodic Laplace problem in three dimensions,

- C++ and Python interface.

Installation
============

The code is available from https://github.com/paulcazeaux/periodic-laplace.git
To download the code, you need to have the Git version control system (http://git-scm.com) installed.

To build and install the module, you need to install the BEM++ library available from https://github.com/bempp/bempp.
Follow the up-to-date instructions at http://www.bempp.org/installation.html to install version 2.0 of the BEM++ library.

Then run the commands:
-	mkdir build; cd build
-	cmake .. && make install


Contact
=======
paul.cazeaux@epfl.ch