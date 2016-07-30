function [U,RAY] = OpenST_LSM3D_TEST(NI, NJ, NK)

H(1) = 1.0 / (NI - 1);
H(2) = 1.0 / (NJ - 1);
H(3) = 1.0 / (NK - 1);

RCV = [.1 .1 .1];
SRC = [.5 .5 .5];

V = ones(NI, NJ, NK);

TSTEP =  max(H(:)) / max(V(:));
EPS = 0.01 * TSTEP;
MAX_ITER = 10;

tic;[U,c,it] = OpenST_LSM3D(V,SRC,H,EPS,MAX_ITER);t = toc;
fprintf('EIKONAL_EX1 test for OpenST_LSM3D MEX\n');
fprintf('converged: %i; iterations: %i; seconds: %.5f\n',c,it,t);
[L1, L2, LINF, UMIN, UMEAN, UMAX] = ex_check(U,SRC,H);
fprintf('L1: %e; L2: %e; LINF: %e; UMIN: %e; UMEAN: %e; UMAX: %e\n', ...
    L1, L2, LINF, UMIN, UMEAN, UMAX);

RAY = OpenST_BRT3D(U,V,H,TSTEP,RCV,SRC);
fprintf('RAY(1): %e %e %e\n',RAY(1,1),RAY(1,2),RAY(1,3));
fprintf('RAY(end): %e %e %e\n',RAY(end,1),RAY(end,2),RAY(end,3));

end

function [L1, L2, LINF, UMIN, UMEAN, UMAX] = ex_check(U,SRC,H)
NI = size(U,1);
NJ = size(U,2);
NK = size(U,3);
NN = NI * NJ * NK;

UEXACT = zeros(NI,NJ,NK);
for i = 0 : (NI - 1)
    for j = 0 : (NJ - 1)
        for k = 0 : (NK - 1)
            di = SRC(1) - i * H(1);
            dj = SRC(2) - j * H(2);
            dk = SRC(3) - k * H(3);
            UEXACT(i + 1,j + 1,k + 1) = sqrt(di^2 + dj ^2 + dk^2);
        end;
    end;
end;

L1 = sum(abs(U(:) - UEXACT(:))) / NN;
L2 = sum(abs(U(:) - UEXACT(:)).^2) / NN;
LINF = max(abs(U(:) - UEXACT(:)));

UMIN = min(U(:));
UMAX = max(U(:));
UMEAN = mean(U(:));
end
