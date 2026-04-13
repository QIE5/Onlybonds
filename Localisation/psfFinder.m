%% Locate 3 different microbubbles, 1 for each region in the frame. Run 
% the code using a selected reference frame. Then, use 
% the centre coordinates (from centroids) and approximate width and height 
% of the bubbles. Continue to re-run code until all 3 bubbles are 
% approximated. These properties are used to form a psfTemplate for each 
% region. Use the psfTemplates to run crossCorrelation.m

% Frame thresholding
function [frameOut] = thresholding (frameIn)
    frameIn = im2gray(frameIn);
    threshold = prctile(frameIn(:), 95);
    frameOut = frameIn;
    frameOut(frameIn <= threshold) = 0;
end

%% PSF info block
% Insert the respective block of code from psfStorage.m. The code block 
% must include: simvid, refFrame1, refFrame2, at least 1 bubble, psf 
% dimensions for each bubble


%% CEUS Stable1
simvid = VideoReader('Phantom Videos/CEUS_Stable1.mp4');
numframes = simvid.NumFrames;
refFrame1 = read(simvid, 1);
refFrame1 = thresholding(refFrame1);
refFrame1 = im2gray(refFrame1);
refFrame2 = read(simvid, 270);
refFrame2 = thresholding(refFrame2);
refFrame3 = read(simvid, 390);
refFrame3 = thresholding(refFrame3);

displayFrame = refFrame1; % Pick the display frame
figure;
imshow(displayFrame);
impixelinfo;
zoom on;
title ('Display PSF regions')
hold on


% Bubble 1
bubble1 = [362,605]; % refFrame3
plot(bubble1(1), bubble1(2), 'b*');
psfWidth1 = 60;
psfHeight1 = 18;
plot ([bubble1(1) - psfWidth1/2 ,bubble1(1) + psfWidth1/2], [bubble1(2), bubble1(2)], 'r-'); 
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + psfHeight1/2,bubble1(2) - psfHeight1/2], 'g-');

% Bubble 2
bubble2 = [662,604]; % refFrame2
plot(bubble2(1), bubble2(2), 'b*');
psfWidth2 = 70;
psfHeight2 = 24;
plot ([bubble2(1) - psfWidth2/2 ,bubble2(1) + psfWidth2/2], [bubble2(2), bubble2(2)], 'r-');%
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + psfHeight2/2,bubble2(2) - psfHeight2/2], 'g-'); 



% Bubble 3
bubble3 = [1206,665]; % Frame 810, refFrame1
plot(bubble3(1), bubble3(2), 'b*');
psfWidth3 = 26;
psfHeight3 = 14;
plot ([bubble3(1) - psfWidth3/2 ,bubble3(1) + psfWidth3/2], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 70
plot ([bubble3(1), bubble3(1)], [ bubble3(2) + psfHeight3/2,bubble3(2) - psfHeight3/2], 'g-'); % height = 14
%% This function determines the psfTemplate based on psf dimensions
function [psfTemplate, box] = findPsfTemplate (frame, psfWidth, psfHeight, bubbleX, bubbleY)
    % patchWidth = psfWidth*1.5;
    % patchHeight = psfHeight*1.5;
    % patchWidth = psfWidth*1.2;
    % patchHeight = psfHeight*1.2;
    patchWidth = psfWidth;
    patchHeight = psfHeight;
    heightDiff = (patchHeight - psfHeight)/2;
    widthDiff = (patchWidth - psfWidth)/2;
    x = bubbleX - (psfWidth)/2 - widthDiff;
    y = bubbleY - (psfHeight)/2 - heightDiff;
    w = patchWidth;
    h = patchHeight;
    box = [x,y, w, h];
    patch = frame(y:y+h-1, x:x+w-1);
    xSum = sum(patch, 1);  % Sum along rows (vertical sum)
    zSum = sum(patch, 2);    % Sum along columns (horizontal sum)

    zSum = zSum / max(zSum);
    xSum = xSum / max(xSum);

    zHalf = find(zSum > 0.5);
    xHalf = find(xSum > 0.5);

    zFwhm = zHalf(end) - zHalf(1) + 1;
    xFwhm = xHalf(end) - xHalf(1) + 1;
    sigmaZ = zFwhm / (2 * sqrt(2 * log(2))); % The relationship between FWHM and standard deviation
    sigmaX = xFwhm / (2 * sqrt(2 * log(2)));

    x1 = -w/2:1:w/2;
    x2 = -h/2:1:h/2;
    [X1,X2] = meshgrid(x1,x2);
    X = [X1(:) X2(:)];
    % normalPsf = exp(-X.^2/(2*sigmaX^2) - Z.^2/(2*sigmaZ^2));
    % normalPsf = normalPsf/max(normalPsf(:));
    mu = [0 0];
    covarMat = [sigmaX^2 0; 0 sigmaZ^2]; % Set cov (X,Z) = 0 as PSF is symmetrical
    normalPsf = mvnpdf(X,mu,covarMat);
    normalPsf = normalPsf / max(normalPsf);
    psfTemplate = reshape(normalPsf,length(x2),length(x1));
end

% Run the function above, entering psf dimensions and bubble centroids
[psfTemplate1, box1] = findPsfTemplate (refFrame1, psfWidth1, psfHeight1, bubble1(1), bubble1(2));
[psfTemplate2, box2] = findPsfTemplate (refFrame2, psfWidth2, psfHeight2, bubble2(1), bubble2(2));
[psfTemplate3, box3] = findPsfTemplate (refFrame1, psfWidth3, psfHeight3, bubble3(1), bubble3(2));

rectangle('Position', box1, 'EdgeColor', 'g', 'LineWidth', 1);
rectangle('Position', box2, 'EdgeColor', 'g', 'LineWidth', 1);
rectangle('Position', box3, 'EdgeColor', 'g', 'LineWidth', 1);
hold off

%% Centroid block
% Use this function to find the centroid locations of all MBs in 1 frame.
% Select an appropriate MB, and approximate its dimensions in the PSF info
% block


% Locate centroids in the reference frame
% refFrame = refFrame1; % Use the line below to pick an appropriate frame


% figure;
% imshow(refFrame)
% impixelinfo;
% title('Find Bubble Centres')








%lsqcurvefit



%% Discretation calculations, ignore for the most part
% ≥0.5λ in the lateral direction , ≥ 0.25λ in the axial direction.
fc = 2.0*1e6;
c = 1540;
lambda = c /fc; % λ = c / f
pixelSize = 3.3879*1e-5;
lambda/4 > pixelSize; % This returns true, so discretation is enough
