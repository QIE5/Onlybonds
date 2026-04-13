function [beamformedImage, srImg, velImg] = ULMPipeline (rawSig, bubbleVid, ...
    constantsParam, beamformParam, svdParam, motionCorrectionParam, localisationParam, ...
    trackingParam, velocityParam)
    % numFrames = 21;
    numFrames = bubbleVid.NumFrames;
    frameRate = constantsParam.frameRate;  % Hz 
    firstFrame = readFrame (bubbleVid);
    firstFrame = im2gray(firstFrame);

    [rows, cols] = size(firstFrame);

    frames = zeros(rows, cols, numFrames, 'like', firstFrame);
    
    for n = 1:numFrames
        frame = read(bubbleVid, n);
        frame = im2gray(frame);
        frames(:, :, n) = frame;
    end
    if strcmp(beamformParam.method, 'DAS')
        beamformedImage = beamform(rawSig, beamformParam, 'SignalType', beamformParam.signal);
    else
        beamformedImage = '';
    end
    if strcmp(svdParam.method, 'SVD')
        % [stack_moving, db_moving] = run_moving_filter(beamformed_data(:,:,1:numFrames));
        % frames = im2uint8(mat2gray(db_moving, [-40 0]));
        stack_filt = svdClutterFilter(frames,4,0);
        stack_filt = svdClutterFilterVideo(frames, 4,0);
    else
        stack_filt = 0;
    end
    % correctedFrames = motionCorrection(bubbleVid);
    if strcmp(motionCorrectionParam.method, 'Motion Correction')
        frames = motionCorrection(frames);
        % numFrames = size(correctedFrames, 4);
    end
    if strcmp(localisationParam.method, 'Cross Correlation')
        localisedBubbleCoords = cell(numFrames, 1);
        for n = 1:numFrames;
            frame = frames(:, :, n);
            % frame = im2gray(frame);
            [localisedBubbleCoords{n}] = crossCorrelation(frame, localisationParam);
        end
    end
    if all(cellfun(@isempty, localisedBubbleCoords))
        error('No bubbles have been localised')
    end
    if strcmp(trackingParam.method, 'Hungarian')
        [tracks, adjacency_tracks, A] = simpletracker(localisedBubbleCoords, ...
        'Method', 'Hungarian');
    elseif strcmp(trackingParam.method, 'Kalman and Hungarian')
        [tracks, adjacency_tracks, A, trackStructs, offsets] = kalmanHungarianTracker1(localisedBubbleCoords, ...
            'MaxGapClosing', 3, ...
            'MaxLinkingDistance', 30, ...
            'ProcessNoise', 10, ...
            'MeasurementNoise', 4, ...
            'InitialVelocityVariance', 100, ...
            'MinTrackLength', 2);
    end
    if strcmp(velocityParam.method, 'Velocity')
        
        [velocity, location, weight] = estimateTrackVelocity(trackStructs, ...
            localisedBubbleCoords, offsets, frameRate);
    end
    [srImg, velImg] = mapping(frame, localisedBubbleCoords, adjacency_tracks, velocity, location, weight, frameRate);
end

%% Pipeline folders
rootDir = fileparts(mfilename('fullpath'));

addpath(fullfile(rootDir, 'Beamforming'));
addpath(fullfile(rootDir, 'Localisation'));
addpath(fullfile(rootDir, 'Mapping'));
addpath(fullfile(rootDir, 'Motion Correction'));
addpath(fullfile(rootDir, 'SVD'));
addpath(fullfile(rootDir, 'Tracking'));

%% ULM constants
constantsParam = struct();
% constantsParam.pixelSize = 3.3879*1e-5; % m % Simulation
constantsParam.pixelSize = 3.0800e-05; % m % Phantom
% constantsParam.frameRate = 400; % Hz, Simulation
constantsParam.frameRate = 40; % Phantom

%% Beamforming parameters: Transducer parameters   
beamformParam = struct();
% beamformParam.method = 'DAS';
beamformParam.method = 'None';
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

%% SVD: parameters
% svdParam.method = 'SVD';
svdParam.method = 'None';




%% Motion correction parameters
% motionCorrectionParam.method = 'Motion Correction';
motionCorrectionParam.method = 'None';

%% Localisation parameters 

if ~exist('psfTemplate1', 'var')
    error('Please extract PSF first, and define microbubbles using psfFinder.m');
end
% localisationParam.method = 'Cross Correlation';
localisationParam.method = 'Cross Correlation';
localisationParam.psfTemplates = {psfTemplate1, psfTemplate2, psfTemplate3};
%% Tracking parameters

% trackingParam.method = 'Hungarian';
trackingParam.method = 'Kalman and Hungarian';
% trackingParam.method = 'None';
%% Velocity parameters
velocityParam.method = 'Velocity';
% velocityParam.method = 'None';
%% Run the pipeline
load("RcvData.mat");
rawSig = RcvData;
% bubbleVid = VideoReader('simulation.mp4');
% bubbleVid = VideoReader('static_background_clutter_filterd.mp4');
bubbleVid = VideoReader('Phantom Videos/CEUS_Stable1.mp4');
[bfImageDB, srImg, velImg] = ULMPipeline (rawSig, bubbleVid, ...
    constantsParam, beamformParam, svdParam, motionCorrectionParam, ...
    localisationParam, trackingParam, velocityParam);




pixelSize = constantsParam.pixelSize;

validMask = velImg.validMask;

%% ============================================================
% 1. PHYSICAL COORDINATE SYSTEM
% ============================================================

W = size(srImg, 2);
H = size(srImg, 1);

dx = pixelSize;
dz = pixelSize;

xAxis_mm = ((1:W) - 0.5) * dx * 1e3;
zAxis_mm = ((1:H) - 0.5) * dz * 1e3;

%% ============================================================
% 2. STRUCTURE IMAGE
% ============================================================

figure;
imagesc(xAxis_mm, zAxis_mm, log(srImg + 1));
axis image off;
colormap(gca, hot);

hold on;

barLength_mm = 10;

xMin_mm = min(xAxis_mm);
xMax_mm = max(xAxis_mm);
zMin_mm = min(zAxis_mm);
zMax_mm = max(zAxis_mm);

x0 = xMin_mm + 0.05 * (xMax_mm - xMin_mm);
y0 = zMax_mm - 0.05 * (zMax_mm - zMin_mm);

plot([x0, x0 + barLength_mm], [y0, y0], 'w', 'LineWidth', 3);
plot([x0, x0], [y0, y0 - barLength_mm], 'w', 'LineWidth', 3);

hold off;

%% ============================================================
% 3. SPEED (mm/s)
% ============================================================

figure;

vx_mm = velImg.vx * dx * 1e3;
vz_mm = velImg.vy * dz * 1e3;

speed = hypot(vx_mm, vz_mm);

maxS = prctile(speed(validMask), 99);
if isempty(maxS) || maxS == 0
    maxS = 1;
end

normalized = speed / maxS;
normalized = max(0, min(1, normalized));

cmap = turbo(256);
idx = round(normalized * 255) + 1;
RGBspeed = ind2rgb(idx, cmap);
RGBspeed(repmat(~validMask, [1 1 3])) = 0;

image(xAxis_mm, zAxis_mm, RGBspeed);
axis image off;
colormap(gca, turbo(256));
set(gca, 'CLim', [0 maxS]);
colorbar;

%% ============================================================
% 4. DIRECTION MAP
% ============================================================

figure;
Hdir = (velImg.direction + pi) / (2*pi);
S = ones(size(Hdir));
V = double(validMask);
RGB = hsv2rgb(cat(3, Hdir, S, V));
image(xAxis_mm, zAxis_mm, RGB);
axis image off;

%% ============================================================
% 5. COLOUR WHEEL
% ============================================================

figure;

nw = 200;
[x,y] = meshgrid(linspace(-1,1,nw));

r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);

Hw = (theta + pi)/(2*pi);
Sw = ones(size(Hw));
Vw = double(r <= 1);

wheel = hsv2rgb(cat(3,Hw,Sw,Vw));

h = image(wheel);
axis image off;
set(h, 'AlphaData', double(r <= 1));


