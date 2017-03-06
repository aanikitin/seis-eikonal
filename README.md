# SEIS-EIKONAL

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-brightgreen.svg)](https://github.com/aanikitin/seis-eikonal/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.269001.svg)](https://doi.org/10.5281/zenodo.269001)

**work in progress**

Eikonal equation solver for calculating first-arrival travel times of seismic 
waves based on the proposed parallel block sweeping methods. These methods are 
based on domain decomposition approach to increase the efficiency of CPU cache 
use. This allows us to achieve high parallel implementation efficiency of 
85-95% (as tested on up to 12 CPU cores). For more information see conference 
abstract [1], extended paper currently submitted for publication.

Recommended build environment is Linux (e.g. Ubuntu 16.04) with a C compiler 
that supports OpenMP 4.0 standard, such as gcc 4.9.0 and higher, and CMake 
3.6.2 and higher. To build the library and test programs, use compile.sh script 
in the project's root directory.

This is a work in progress, the library will be updated in the future to 
improve its interface and performance and to add documentation and new 
functionality. Latest versions are available at 
https://github.com/aanikitin/seis-eikonal.

Authors:<br />
Alexandr Nikitin -  PhD student, junior researcher, IPGG SB RAS - parallel 
algorithm and software development.<br />
Anton Duchkov, PhD - head of laboratory, IPGG SB RAS - PhD thesis advisor 
(geophysics).<br />
Alexandr Serdyukov, PhD - researcher, IPGG SB RAS - consultant (numerical 
methods).<br />

References:<br />
1) Nikitin Alexandr A., Serdyukov Alexandr S., Duchkov Anton A. <a 
href="http://elibrary.ru/item.asp?id=25994951">Optimization of parallel 
sweeping methods of numerical computation of seismic wave travel times for 
shared memory computing systems</a> // ИНТЕРЭКСПО 
ГЕО-СИБИРЬ. – 2016. – V. 2. – N. 1. – P. 241-245.<br />

Last change: 2017-03-07
