data = load("RcvData.mat");
rf = data.RcvData;

%% Transducer parameters
param.fc = 2.841e6;           % probe central frequency (Hz)
param.fs = 14.205e6;          % sampling rate (Hz)
param.c = 1540;               % speed of sound (m/s)
param.pitch = 2.7e-4;         % element pitch (m)
param.fnumber = 1.4;          % f-number used in paper

%% Image reconstruction parameters
pixelmap_x = -0.06:(2.7*1e-4):0.06; 
pixelmap_z = 0:(3.4*1e-4):0.08;     

[X, Z] = meshgrid(pixelmap_x, pixelmap_z);

%% Virtual source
vsource = [0, -10000]; % x0=steering_angle=0, z0 is large and negative 

%% Beamforming

[bfSIG, M] = ezdas(rf, X(:), Z(:), vsource, param);

% Reshape to image dimensions
bfImage = reshape(bfSIG, size(X));

%% Display the result
figure('Position', [100, 100, 800, 600]);

% Display in dB scale
bfImage_dB = 20*log10(abs(bfImage) + eps);
bfImage_dB = bfImage_dB - max(bfImage_dB(:));  % Normalize to 0 dB

imagesc(pixelmap_x*1e3, pixelmap_z*1e3, bfImage_dB, [-60, 0]);
colormap('gray');
colorbar;
xlabel('Lateral Distance (mm)');
ylabel('Axial Distance (mm)');
title('DAS Beamformed Image');
axis equal tight;

fprintf('Beamforming complete!\n');
