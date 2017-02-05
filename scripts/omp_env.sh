#!/bin/bash
# script used to setup OpenMP environment for testing

export OMP_DISPLAY_ENV=verbose
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
