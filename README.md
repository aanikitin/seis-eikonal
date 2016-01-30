# SEIS-EIKONAL 0.1-POC

Proof of Concept version

Eikonal equation solver for calculating first-arrival traveltimes of seismic 
waves currently based on parallel implementation of the FSM algorithm. Two new 
implementations are proposed for shared memory architectures based on 3D grid 
decomposition into tasks to achieve high performance. One based on static 
assignment of tasks to threads using manual scheduling. The other is using task 
depend clause introduced in OpenMP 4.0 to avoid the need for global 
synchronization of task execution. Preliminary testing demonstrates that both 
approaches show similar performance under light system load, with the later 
approach becoming faster during increased system load.

This is a work in progress. We are currently preparing an article with detailed 
description of the proposed algorithms and in-depth testing, and are planning 
to submit it for review and publication in the first half of 2016. In the 
meantime new major versions of the project will be uploaded to Zenodo.org 
repository. **If you are planning to use our proposed algorithms in your 
research before we had a chance to publish an article, please cite our research 
using DOIs for our source code uploaded to Zenodo**. Thank you.

Authors:
Alexandr Nikitin -  PhD student, junior researcher, IPGG SB RAS - parallel 
algorithm and software development. 
Anton Duchkov, PhD - head of laboratory, IPGG SB RAS - PhD thesis advisor 
(geophysics). 
Alexandr Serdyukov, PhD - researcher, IPGG SB RAS - consultant (numerical 
methods).
