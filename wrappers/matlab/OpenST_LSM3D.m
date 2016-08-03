function [U,varargout] = OpenST_LSM3D(V,SRC,varargin)

if nargin > 2
    H = varargin{1};
else
    H = [1, 1, 1];
end;

if nargin > 3
    EPS = varargin{2};
else
    EPS = 0.01 * min(H(:)) / max(V(:));
end;

if nargin > 4
    MAX_ITER = varargin{3};
else
    MAX_ITER = 100;
end;

V = permute(V, [3 2 1]);
[U,c,it] = OpenST_LSM3D_mex(V,SRC,H,EPS,MAX_ITER);
U = permute(U, [3 2 1]);

if nargout > 1
    varargout{1} = c;
end;

if nargout > 2
    varargout{2} = it;
end;

end
