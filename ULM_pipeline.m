%%

function [beamformedImage, SR_img, BW] = ULMPipeline (rawSig, bubbleVid, beamformParam, localisationParam)
    numFrames = 10;
    beamformedImage = beamform(rawSig, beamformParam, 'SignalType', beamformParam.signal);
    correctedFrames = motionCorrection(bubbleVid);
    % numFrames = size(correctedFrames, 4);
    localisedBubbleCoords = cell(numFrames, 1);
    for n = 1:numFrames;
        frame = correctedFrames(:, :, :, n);
        [localisedBubbleCoords{n}, boxes] = localisation(frame, localisationParam);
    end
    [tracks, adjacency_tracks, A] = simpletracker(localisedBubbleCoords, ...
    'Method', 'Hungarian');
    [SR_img, BW] = mapping(frame, localisedBubbleCoords, tracks, adjacency_tracks, A);
end


%% Pipeline folders
rootDir = fileparts(mfilename('fullpath'));

addpath(fullfile(rootDir, 'Beamforming'));
addpath(fullfile(rootDir, 'Localisation'));
addpath(fullfile(rootDir, 'Mapping'));
addpath(fullfile(rootDir, 'Motion Correction'));
addpath(fullfile(rootDir, 'Tracking'));


%% Beamforming parameters: Transducer parameters   
beamformParam.type = 'das';
beamformParam.signal = 'RF';
beamformParam.fc = 2.0*1e6;           % transmission frequency
beamformParam.fs = 14.205e6;          % sampling rate (Hz)
beamformParam.c = 1540;               % speed of sound (m/s) - this can be calculated using the speed of sound function in the paper
                              % the paper does make the argument that
                              % optimising the speed of sound may be
                              % valuable in SRUS
beamformParam.pitch = 2.7e-4;         % element pitch (m)
% transducerParam.fnumber = 1.4;          % f-number used in paper

% f-number calculation
% elementWidth = 2.3*1e-4; 
% lambda = transducerParam.c / transducerParam.fc;  % Calculate wavelength based on speed of sound and central frequency
% f = @(th,width,lambda)... 
%     abs(cos(th)*sinc(width/lambda*sin(th))-0.71); 
% alpha = fminbnd(@(th) f(th,elementWidth,lambda),0,pi/2); 
% transducerParam.fnumber = 1/2/tan(alpha); % f-number = 0.5786

beamformParam.fnumber = 0.5786; % for simplification, users should enter their f-number

%% Beamforming parameters: Image reconstruction parameters

beamformParam.pixelmapX = -0.06:(2.7*1e-4):0.06; 
beamformParam.pixelmapZ = 0:(3.4*1e-4):0.08;     

%% Localisation parameters 
localisationParam.threshold = 99.5;
localisationParam.psfHeight = 33;
localisationParam.psfWidth = 183;
%%

load("RcvData.mat");
rawSig = RcvData;
bubbleVid = VideoReader('simulation.mp4');
[bfImageDB, SR_img, BW] = ULMPipeline (rawSig, bubbleVid, beamformParam, localisationParam);

%% Display the result
figure('Position', [100, 100, 800, 600]);

% Display in dB scale

imagesc(beamformParam.pixelmapX, beamformParam.pixelmapZ, bfImageDB, [-40, 0]); 
colormap('gray');
colorbar;
xlabel('Lateral Distance (m)');
ylabel('Axial Distance (m)');
title('DAS Beamformed Image');
axis image;
%% Final image
figure;
subplot(1,2,1);
imagesc(log(SR_img + 1));
axis image off; colormap hot;
title('ULM Density Map');

subplot(1,2,2);
imshow(BW);
title('Reconstructed Vessel');