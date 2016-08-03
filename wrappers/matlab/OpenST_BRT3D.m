function RAY = OpenST_BRT3D(T,V,H,TSTEP,RCV,SRC,varargin)
% RAY = OpenST_BRT3D(T,V,H,TSTEP,RCV,SRC);

if nargin > 6
    MAX_SEG = varargin{1};
else
    MAX_SEG = numel(T);
end;

T = permute(T, [3 2 1]);
V = permute(V, [3 2 1]);
RAY = OpenST_BRT3D_mex(T,V,H,TSTEP,RCV,SRC,MAX_SEG);
RAY = permute(RAY, [2 1]);

end
