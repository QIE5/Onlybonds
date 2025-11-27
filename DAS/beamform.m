data = load("RcvData.mat");
rf = data.RcvData;

%% Transducer parameters
param.fc = 2.841e6;           % probe central frequency (Hz)
param.fs = 14.205e6;          % sampling rate (Hz)
param.c = 1540;               % speed of sound (m/s) - this can be calculated using the speed of sound function in the paper
                              % the paper does make the argument that
                              % optimising the speed of sound may be
                              % valuable in SRUS
param.pitch = 2.7e-4;         % element pitch (m)
% param.fnumber = 1.4;          % f-number used in paper

% instead, I will use the matlab script from the paper to calculate the
% f-number
steering_angle = 0;
element_width = 2.3*1e-4; 
lambda = param.c / param.fc;  % Calculate wavelength based on speed of sound and central frequency
f = @(th,width,lambda)... 
    abs(cos(th)*sinc(width/lambda*sin(th))-0.71); 
alpha = fminbnd(@(steering_angle) f(steering_angle,element_width,lambda),0,pi/2); 
param.fnumber = 1/2/tan(alpha) % f-number = 0.6532

%% Image reconstruction parameters
pixelmap_x = -0.06:(2.7*1e-4):0.06; 
pixelmap_z = 0:(3.4*1e-4):0.08;     

[X, Z] = meshgrid(pixelmap_x, pixelmap_z);

%% Virtual source
vsource = [0, -10000]; % x0=steering_angle=0, z0 is large and negative 
% if there is no steering angle, then the virtual source is directly above
% the array. This makes theta perpendicular to the z axis
% beta can also be calculated based on the virtual source

%% Beamforming

[bfSIG, M] = ezdas(rf, X(:), Z(:), vsource, param);
% This function doesn't apply directivity weighting. Directivity = 1 in our
% case due to the steering angle of 0.

% Reshape to image dimensions
bfImage = reshape(bfSIG, size(X));

%% Display the result
figure('Position', [100, 100, 800, 600]);
% Display in dB scale
bfImage_dB = 20*log10(abs(bfImage) + eps); % b-mode is obtained via logging. Image looks strange without it
                                           %abs(bfImage) takes the
                                           %magnitude
bfImage_dB = bfImage_dB - max(bfImage_dB(:));  % Normalise to 0 dB

imagesc(pixelmap_x, pixelmap_z, bfImage_dB, [-60, 0]); % Axes are in mm
colormap('gray');
colorbar;
xlabel('Lateral Distance (m)');
ylabel('Axial Distance (m)');
title('DAS Beamformed Image');
axis image;

