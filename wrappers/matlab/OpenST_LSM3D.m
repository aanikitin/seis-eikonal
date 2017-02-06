function [U,varargout] = OpenST_LSM3D(V,SRC,varargin)
%OpenST_LSM3D Solve the eikonal equation on a cartesian grid.
%   U = OpenST_LSM3D(V,SRC) solves the eikonal equation for the 
%   three-dimensional velocity model V(i,j,k) and the source located at 
%   coordinates SRC = [ci,cj,ck]. U is the unknown such as first-arrival 
%   travel times in seismic applications. Grid steps default to H = 
%   [1,1,1]. The origin of the cartesian coordinate system corresponds to 
%   grid point [i,j,k] = [0,0,0]. Coordinates of all grid points are 
%   calculated as [(i - 1)*H(1),(j - 1)*H(2),(k - 1)*H(3)].
%
%   U = OpenST_LSM3D(V,SRC,H) will set grid steps to H = [HI,HJ,HK].
%
%   U = OpenST_LSM3D(V,SRC,H,EPS) will set the stopping criterion parameter 
%   for the locking sweeping numerical method to EPS. The iterative method 
%   stops when the absolute maximum between solutions obtained at the 
%   current and previous iterations becomes less than EPS (L-infinity 
%   norm). Defaults to the ratio of minimum grid step to maximum velocity.
%
%   U = OpenST_LSM3D(V,SRC,H,EPS,MAX_ITER) will set the maximum number of 
%   iterations to MAX_ITER. Defaults to 100.
%
%   [U,C] = OpenST_LSM3D(...) returns boolean variable C that specifies
%   whether solution converged to EPS.
%
%   [U,C,IT] = OpenST_LSM3D(...) returns the number of iterations 
%   completed in IT.
%
%   Visit https://github.com/aanikitin/seis-eikonal for latest version.
%
%   See also OpenST_MEX_Compile, OpenST_BRT3D.

%   Copyright 2014-2017 Alexandr Nikitin.

if nargin > 2
    H = varargin{1};
else
    H = [1, 1, 1];
end;

if nargin > 3
    EPS = varargin{2};
else
    EPS = 1.0 * min(H(:)) / max(V(:));
end;

if nargin > 4
    MAX_ITER = varargin{3};
else
    MAX_ITER = 100;
end;

V = permute(V, [3 2 1]);
[U,C,IT] = OpenST_LSM3D_mex(V,SRC,H,EPS,MAX_ITER);
U = permute(U, [3 2 1]);

if nargout > 1
    varargout{1} = C;
end;

if nargout > 2
    varargout{2} = IT;
end;

end
