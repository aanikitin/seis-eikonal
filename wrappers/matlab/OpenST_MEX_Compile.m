function OpenST_MEX_Compile(varargin)
%OpenST_MEX_Compile Compiles OpenST MEX functions.
%
%   Visit https://github.com/aanikitin/seis-eikonal for latest version.
%
%   See also OpenST_LSM3D, OpenST_BRT3D.

%   Copyright 2014-2017 Alexandr Nikitin.

if nargin > 0
    VERBOSE = varargin{1};
else
    VERBOSE = false;
end;

% path to OpenST static library
OpenST_static_lib_path = '../../lib/Release';
% fullpath to OpenST static library
OpenST_static_lib_fullpath = sprintf('%s/openst-1-static.lib', ...
    OpenST_static_lib_path);

if VERBOSE
    mex('-v','-largeArrayDims',...
        'COMPFLAGS=$COMPFLAGS /I ../../include',...
        'OpenST_LSM3D_mex.c','OpenST_MEX.c',OpenST_static_lib_fullpath);
    mex('-v','-largeArrayDims',...
        'COMPFLAGS=$COMPFLAGS /I ../../include',...
        'OpenST_BRT3D_mex.c','OpenST_MEX.c',OpenST_static_lib_fullpath);
else
    mex('-largeArrayDims','COMPFLAGS=$COMPFLAGS /I ../../include',...
        'OpenST_LSM3D_mex.c','OpenST_MEX.c',OpenST_static_lib_fullpath);
    mex('-largeArrayDims','COMPFLAGS=$COMPFLAGS /I ../../include',...
        'OpenST_BRT3D_mex.c','OpenST_MEX.c',OpenST_static_lib_fullpath);
end;

end
