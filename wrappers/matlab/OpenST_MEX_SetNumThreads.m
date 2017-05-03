function NumThreadsSet = OpenST_MEX_SetNumThreads(NumThreads)
clear OpenST_LSM3D_mex OpenST_BRT3D_mex;
setenv('OMP_NUM_THREADS',sprintf('%i',NumThreads));
NumThreadsSet = getenv('OMP_NUM_THREADS');
end
