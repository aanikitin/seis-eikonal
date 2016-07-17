# SEIS-EIKONAL

[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-brightgreen.svg)](https://github.com/aanikitin/seis-eikonal/blob/master/LICENSE)

**work in progress, not for production use**

Eikonal equation solver for calculating first-arrival traveltimes of seismic 
waves currently based on parallel implementation of the FSM algorithm. Two new 
implementations are proposed for shared memory architectures based on 3D grid 
decomposition into tasks to achieve high performance. One based on static 
assignment of tasks to threads using manual scheduling. The other is using task 
depend clause introduced in OpenMP 4.0 to avoid the need for global 
synchronization of task execution. Preliminary testing demonstrates that both 
approaches show similar performance under light system load, with the later 
approach becoming faster during increased system load.

We are currently preparing an article with detailed 
description of the proposed algorithms and in-depth testing, and are planning 
to submit it for review and publication in the first half of 2016.

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
