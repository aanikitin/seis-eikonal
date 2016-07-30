function RAY = OpenST_BRT3D(T,V,H,TSTEP,RCV,SRC)
% RAY = OpenST_BRT3D(T,V,H,TSTEP,RCV,SRC);

T = permute(T, [3 2 1]);
V = permute(V, [3 2 1]);
RAY = OpenST_BRT3D_mex(T,V,H,TSTEP,RCV,SRC);
RAY = permute(RAY, [2 1]);

end
