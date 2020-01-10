# Separation of Shape and Temporal Patterns in Time Series

This repository contains a Python wrapper for the C++ code for averaging time series using elastic kernels defined and used in [1] with some extensions dedicated to the separation of shape and temporal pattern in time series as described in [2]. 
This repository is self contained and extends [eKATS](https://github.com/pfmarteau/eKATS) repository.

## Installation (tested on Ubuntu 16.04, 18.04)

To compile and install the C++ code with its python wrapper:

$ sh install.sh

## Running the noisy ellipses example

$ Python3 STS2_noisyEllipses.py


## References

[1] Marteau, P.F., Times series averaging and denoising from a probabilistic perspective on time-elastic kernels International Journal of Applied Mathematics and Computer Science, De Gruyter [pdf](https://arxiv.org/abs/1611.09194)

[2] Marteau, P.F., On the separation of shape and temporal patterns in time series -Application to signature authentication-  [pdf](https://arxiv.org/abs/1911.09360)
