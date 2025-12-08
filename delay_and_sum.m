function [bfSIG, M] = ezdas(SIG, x, z, vsource, param)
% EZDAS Simple Delay-and-Sum beamformer
%
% INPUTS:
%   SIG     - RF or IQ data [nSamples x nElements]
%   x       - x-coordinates of image pixels (vector or grid)
%   z       - z-coordinates of image pixels (vector or grid)
%   vsource - [x0, z0] source position
%   param   - struct with fields:
%               pitch     : element pitch (m)
%               c         : speed of sound (m/s)
%               fs        : sampling frequency (Hz)
%               fnumber   : F-number for aperture
%               fc        : center frequency (Hz)
%
% OUTPUTS:
%   bfSIG   - beamformed image [length(x) x length(z)]
%   M       - sparse interpolation matrix

% Ensure column vectors
x = x(:);
z = z(:);
siz0 = size(x);

[n1, nc] = size(SIG);  % nSamples x nElements

% Element positions
xe = ((0:nc-1) - (nc-1)/2) * param.pitch;
L = xe(end) - xe(1);

x0 = vsource(1); z0 = vsource(2);

% Distances
dTX = hypot(x - x0, z - z0) - hypot((abs(x0)-L/2)*(abs(x0)>L/2), z0);
dRX = hypot(x - xe', z);    % xe' makes broadcasting work

% Time-of-flight (seconds)
tau = (dTX + dRX)/param.c;

% Sample indices
idxt = tau * param.fs + 1;

% Valid indices
I = idxt >= 1 & idxt <= n1-1;

% Aperture mask
Iaperture = abs(x - xe') <= (z/2/param.fnumber);
I = I & Iaperture;

% Flatten indices for interpolation
idx = idxt + (0:nc-1)*n1;        % linear indexing for each element
idx = idx(I);
idxf = floor(idx);
frac = idx - idxf;                % fractional delay

% Row indices for sparse matrix
[i, ~] = find(I);

% Linear interpolation weights
s = [1 - frac; frac];

% Phase compensation for complex IQ
if ~isreal(SIG)
    s = s .* exp(2i*pi*param.fc * [tau(I); tau(I)]);
end

% Build sparse interpolation matrix
M = sparse([i; i], [idxf; idxf+1], s, numel(x), n1*nc);

% Beamformed image
bfSIG = reshape(M * SIG(:), siz0);

end
