function [bfSIG,M] = ezdas(SIG,x,z,vsource,param)

siz0 = size(x);
[n1,nc] = size(SIG);
x = x(:); z = z(:);

xe = ((0:nc-1)-(nc-1)/2)*param.pitch;
L = xe(end)-xe(1);

x0 = vsource(1); z0 = vsource(2);

dTX = hypot(x-x0, z-z0) - ...
    hypot((abs(x0)-L/2) * (abs(x0)>L/2),z0);
dRX = hypot(x-xe,z);

tau = (dTX+dRX)/param.c;

idxt = tau*param.fs + 1;

I = idxt>=1 & idxt<=n1-1;
Iaperture = abs(x-xe)<=(z/2/param.fnumber);
I = I&Iaperture;

idx = idxt + (0:nc-1)*n1;
idx = idx(I);
idxf = floor(idx);
idx = idxf-idx;

[i,~] = find(I);
s = [idx+1;-idx];
if ~isreal(SIG)
    s = s .* exp(2i*pi*param.fc*[tau(I); tau(I)]);
end
M = sparse([i;i], [idxf;idxf+1],s,numel(x),n1*nc);

bfSIG = reshape(M*SIG(:), siz0);