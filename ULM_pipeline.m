function [beamformedImage, srImg, velImg] = ULMPipeline (rawSig, bubbleVid, ...
    beamformParam, svdParam, motionCorrectionParam, localisationParam, ...
    trackingParam, velocityParam)
    numFrames = 21;
    numFrames = bubbleVid.NumFrames;
    frameRate = bubbleVid.FrameRate;
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


%% Beamforming parameters: Transducer parameters   
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
bubbleVid = VideoReader('static_background_clutter_filterd.mp4');
[bfImageDB, srImg, velImg] = ULMPipeline (rawSig, bubbleVid, beamformParam, ...
    svdParam, motionCorrectionParam, localisationParam, trackingParam, ...
    velocityParam);

% %% Display the result
% figure('Position', [100, 100, 800, 600]);
% 
% % Display in dB scale
% 
% image = imagesc(beamformParam.pixelmapX, beamformParam.pixelmapZ, bfImageDB, [-40, 0]); 
% colormap('gray');
% colorbar;
% xlabel('Lateral Distance (m)');
% ylabel('Axial Distance (m)');
% title('DAS Beamformed Image');
% axis image;
%% Final image
figure;
validMask = velImg.validMask;

%% Structure
subplot(2,3,1);
imagesc(log(srImg + 1));
axis image off; colormap(gca, hot);
title('Density');

hold on;

% --- Scale bar (20 pixels = 1 µm) ---
barWidthPx = (1e-3) / pixelSize;
% barLength = 20;  % pixels

% Position (bottom-left corner)
xStart = 15;
yStart = size(srImg,1) - 10;

% Draw line
plot([xStart, xStart + barWidthPx], [yStart, yStart], ...
    'w', 'LineWidth', 2);

% % Label
% text(xStart + barLength/2, yStart - 40, '1 \mum', ...
%     'Color', 'w', ...
%     'HorizontalAlignment', 'center', ...
%     'VerticalAlignment', 'top', ...
%     'FontSize', 8);

hold off;
%% Speed
subplot(2,3,2);
maxS = prctile(velImg.speed(validMask), 99);
if isempty(maxS) || maxS == 0, maxS = 1; end
h = imagesc(velImg.speed, [0 maxS]);
axis image off; 
colormap(gca, turbo(256));
set(h, 'AlphaData', double(validMask));
colorbar; title('Speed');

%% Velocity X
subplot(2,3,3);
maxVx = prctile(abs(velImg.vx(validMask)), 99);
if isempty(maxVx) || maxVx == 0, maxVx = 1; end
h = imagesc(velImg.vx, [-maxVx maxVx]);
axis image off; set(gca, 'Color', 'k');
colormap(gca, parula(256));
set(h, 'AlphaData', double(validMask));
colorbar; title('Velocity X');

%% Velocity Y
subplot(2,3,4);
maxVy = prctile(abs(velImg.vy(validMask)), 99);
if isempty(maxVy) || maxVy == 0, maxVy = 1; end
h = imagesc(velImg.vy, [-maxVy maxVy]);
axis image off; set(gca, 'Color', 'k');
colormap(gca, parula(256));
set(h, 'AlphaData', double(validMask));
colorbar; title('Velocity Y');

%% Direction
subplot(2,3,5);
H = (velImg.direction + pi) / (2*pi);
S = ones(size(H));
V = double(validMask);
RGB = hsv2rgb(cat(3, H, S, V));
image(RGB); axis image off;
title('Direction');

%% Colour wheel
subplot(2,3,6);
nw = 200;
[x,y] = meshgrid(linspace(-1,1,nw));
r = sqrt(x.^2 + y.^2);
theta = atan2(y,x);
Hw = (theta+pi)/(2*pi);
Sw = ones(size(Hw));
Vw = double(r <= 1);
wheel = hsv2rgb(cat(3,Hw,Sw,Vw));
image(wheel); axis image off;
title('Direction key', 'FontSize', 8);


