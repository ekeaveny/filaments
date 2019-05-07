# 3D Matlab RPY code: README

Supplementary Matlab/Octave code to ['Methods for suspensions of passive and active filaments'](https://arxiv.org/abs/1903.12609), 2019, by SF Schoeller, AK Townsend, TA Westwood & EE Keaveny.

## 1. Whom do I talk to?
* For contact details, see README.md in the [parent directory](https://github.com/ekeaveny/filaments/).

## 2. What does this code do?
This code demonstrates the use of the method described in [the paper](https://arxiv.org/abs/1903.12609) in simulating a single flexible filament falling under gravity in an infinite domain.

It uses the 'EJBb' version of Broyden's method (Algorithm 2 in the paper) with a reduced 'robot arm' system of nonlinear equations. For the hydrodynamic solver, it uses the RPY tensors.

## 3. How do I run it?
To use, just run `main.m` in Matlab.

## 4. How should I best understand it?
The code is fully commented and docstringed, but we recommend understanding the simpler 2D code first.

## 5. How does it differ from the 2D code?
Structurally, the 3D code is very similar to the 2D code, with one important exception. Unlike the 2D code, the 3D code is object oriented, with many functions and properties existing inside the Filament object. However, you should find that the order of the code is fundamentally the same between the 2D and 3D versions, and care has been taken so that equivalent functions and properties have the same names, and the same comments, in each version.
