function [RAY,varargout] = OpenST_BRT3D(T,V,H,TSTEP,RCV,SRC,varargin)
%OpenST_BRT3D Perform back ray tracing for the source-receiver pair.
%
%   Visit https://github.com/aanikitin/seis-eikonal for latest version.
%
%   See also OpenST_MEX_Compile, OpenST_LSM3D.

%   Copyright 2014-2017 Alexandr Nikitin.

if nargin > 6
    MAX_SEG = varargin{1};
else
    MAX_SEG = numel(T);
end;

T = permute(T, [3 2 1]);
V = permute(V, [3 2 1]);
[RAY,COMPTIME] = OpenST_BRT3D_mex(T,V,H,TSTEP,RCV,SRC,MAX_SEG);
RAY = permute(RAY, [2 1]);

if nargout > 1
    varargout{1} = COMPTIME;
end;

end
