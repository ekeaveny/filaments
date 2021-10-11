# filaments: README

Supplementary Matlab/Octave code to ['Methods for suspensions of passive and active filaments'](https://arxiv.org/abs/1903.12609), 2019, by SF Schoeller, AK Townsend, TA Westwood & EE Keaveny.

## 1. Whom do I talk to?
The authors can be contacted at:

* Simon Schoeller: `simon.schoeller14@imperial.ac.uk`
* Adam Townsend: `adam.townsend@imperial.ac.uk`
* Tim Westwood: `t.westwood16@imperial.ac.uk`
* Eric Keaveny: `e.keaveny@imperial.ac.uk`

## 2. What code is here?
1. [`/matlab-2d-rpy/`](https://github.com/ekeaveny/filaments/tree/master/matlab-2d-rpy) A 2D Matlab/Octave code which uses the method described in [the paper](https://arxiv.org/abs/1903.12609) to simulate a single flexible filament falling under gravity in an unbounded domain.
2. [`/matlab-3d-rpy/`](https://github.com/ekeaveny/filaments/tree/master/matlab-3d-rpy) A 3D Matlab/Octave code which does the same, using quaternions and geometric time integration as described in [the paper](https://arxiv.org/abs/1903.12609) to handle the 3D motion.

## Acknowledgements:
11-OCT-2021: We thank Robert (Bob) Swallow of Trinity University for alerting us to errors in the approximate Jacobian.
