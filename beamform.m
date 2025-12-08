% This uses ezdas.m

data = load("RcvData.mat");
rf = data.RcvData;



%% Transducer parameters   
param.fc = 2.0*1e6;           % transmission frequency
param.fs = 14.205e6;          % sampling rate (Hz)
param.c = 1540;               % speed of sound (m/s) - this can be calculated using the speed of sound function in the paper
                              % the paper does make the argument that
                              % optimising the speed of sound may be
                              % valuable in SRUS
param.pitch = 2.7e-4;         % element pitch (m)
% param.fnumber = 1.4;          % f-number used in paper

% instead, I will use the matlab script from the paper to calculate the
% f-number
elementWidth = 2.3*1e-4; 
lambda = param.c / param.fc;  % Calculate wavelength based on speed of sound and central frequency
f = @(th,width,lambda)... 
    abs(cos(th)*sinc(width/lambda*sin(th))-0.71); 
alpha = fminbnd(@(th) f(th,elementWidth,lambda),0,pi/2); 
param.fnumber = 1/2/tan(alpha) % f-number = 0.5786


%% Image reconstruction parameters
pixelmapX = -0.06:(2.7*1e-4):0.06; 
pixelmapZ = 0:(3.4*1e-4):0.08;     

[X, Z] = meshgrid(pixelmapX, pixelmapZ);

%% Virtual source
vsource = [0, -10000]; % x0=steering_angle=0, z0 is large and negative 
% if there is no steering angle, then the virtual source is directly above
% the array. This makes theta perpendicular to the z axis
% beta can also be calculated based on the virtual source

%% Beamforming
[bfSIG, M] = ezdas(rf, X(:), Z(:), vsource, param);

% Reshape to image dimensions
bfImage = reshape(bfSIG, size(X));
%% Display the result
figure('Position', [100, 100, 800, 600]);
% Display in dB scale
%envelope = abs((bfImage)); %This is the same as the line below:
envelope = hilbert(bfImage);
envelope = abs(envelope);
bfImageDB = 20*log10(envelope / max(envelope(:)) + eps);

imagesc(pixelmapX, pixelmapZ, bfImageDB, [-40, 0]); 
colormap('gray');
colorbar;
xlabel('Lateral Distance (m)');
ylabel('Axial Distance (m)');
title('DAS Beamformed Image');
axis image;

