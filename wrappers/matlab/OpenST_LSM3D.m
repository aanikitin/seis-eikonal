function [U,c,it] = OpenST_LSM3D(V,SRC,H,EPS,MAX_ITER)

V = permute(V, [3 2 1]);
[U,c,it] = OpenST_LSM3D_mex(V,SRC,H,EPS,MAX_ITER);
U = permute(U, [3 2 1]);

end
