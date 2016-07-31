function [U,RAY] = OpenST_LSM3D_TEST(NI, NJ, NK)

% grid step
H(1) = 1.0 / (NI - 1);
H(2) = 1.0 / (NJ - 1);
H(3) = 1.0 / (NK - 1);
% receiver coordinates
RCVK = 0.0;
[RCVJ, RCVI] = meshgrid(.1:.1:.9, .1:.1:.9);
RCV = [RCVI(:), RCVJ(:), RCVK * ones(numel(RCVJ),1)];
% source coordinates
SRC = [.5 .5 .5];
% wave speed
V = ones(NI, NJ, NK);
% back ray tracing time step
TSTEP =  max(H(:)) / max(V(:));
% eikonal stopping criterion
EPS = 0.01 * TSTEP;
% maximum number of iterations for LSM
MAX_ITER = 10;

% compute traveltimes using LSM
tic;[U,c,it] = OpenST_LSM3D(V,SRC,H,EPS,MAX_ITER);t = toc;

fprintf('EIKONAL_EX1 test for OpenST_LSM3D MEX\n');
fprintf('converged: %i; iterations: %i; seconds: %.5f\n',c,it,t);
[L1, L2, LINF, UMIN, UMEAN, UMAX] = ex_check(U,SRC,H);
fprintf('L1: %e; L2: %e; LINF: %e; UMIN: %e; UMEAN: %e; UMAX: %e\n', ...
    L1, L2, LINF, UMIN, UMEAN, UMAX);

hf = figure();
ha = axes('Parent',hf);
xi = 2;
yi = 1;
zi = 3;
x = 0 : H(2) : ((NJ - 1) * H(xi));
y = 0 : H(1) : ((NI - 1) * H(yi));
z = 0 : H(3) : ((NK - 1) * H(zi));
[X,Y,Z] = meshgrid(x,y,z);
hs = slice(X,Y,Z,V,.5,.5,.5,'Parent',ha);
set(hs,'LineStyle','none','FaceAlpha',0.4);
hold(ha,'on');

plot3(ha,SRC(1,xi),SRC(1,yi),SRC(1,zi),'*');
for i = 1:size(RCV,1)
    plot3(ha,RCV(i,xi),RCV(i,yi),RCV(i,zi),'*');
end;

RAY = cell(size(SRC,1) * size(RCV,1),1);
t = zeros(size(SRC,1) * size(RCV,1),1);

% perform back ray tracing
for i = 1:size(RCV,1)
    tic;RAY{i} = OpenST_BRT3D(U,V,H,TSTEP,RCV(i,:),SRC);t(i) = toc;
    plot3(ha,RAY{i}(:,xi),RAY{i}(:,yi),RAY{i}(:,zi),'-');    
end;
fprintf('BRT3D: min %e, mean %e, max %e sec.\n',min(t),mean(t),max(t));
set(ha,'YDir','reverse');
set(ha,'ZDir','reverse');
axis(ha,'square');
xlabel('j');
ylabel('i');
zlabel('k');
hold(ha,'off');
colorbar(ha);

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
