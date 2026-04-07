%% Locate 3 different microbubbles, 1 for each region in the frame. Run 
% the code using a selected reference frame. Then, use 
% the centre coordinates (from centroids) and approximate width and height 
% of the bubbles. Continue to re-run code until all 3 bubbles are 
% approximated. These properties are used to form a psfTemplate for each 
% region. Use the psfTemplates to run crossCorrelation.m

%% PSF info block
% Insert the respective block of code from psfStorage.m. The code block 
% must include: simvid, refFrame1, refFrame2, at least 1 bubble, psf 
% dimensions for each bubble

% Clutter filter video
simvid = VideoReader('static_background_clutter_filterd.mp4');
refFrame1 = read(simvid, 70);
refFrame2 = read(simvid, 140);
figure;
imshow(refFrame1);
title ('Display PSF regions')
hold on

% Bubble 1
bubble1 = [113,44]; % Frame 70
plot(bubble1(1), bubble1(2), 'b*');
plot ([bubble1(1) - 35 ,bubble1(1) + 35], [bubble1(2), bubble1(2)], 'r-'); % PSF width = 70
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + 6,bubble1(2) - 6], 'g-'); % height = 12
psfWidth1 = 70;
psfHeight1 = 12;

% Bubble 2
bubble2 = [94,151]; % Frame 70
plot(bubble2(1), bubble2(2), 'b*');
% Bubble 2, frame 140
plot ([bubble2(1) - 34 ,bubble2(1) + 34], [bubble2(2), bubble2(2)], 'r-'); % PSF width = 64
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 6,bubble2(2) - 6], 'g-'); % height = 12
psfWidth2 = 64;
psfHeight2 = 12;

% Bubble 3
bubble3 = [82,244]; % Frame 70
plot(bubble3(1), bubble3(2), 'b*');
plot ([bubble3(1) - 34 ,bubble3(1) + 34], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 68
plot ([bubble3(1), bubble3(1)], [ bubble3(2) + 6,bubble3(2) - 5], 'g-'); % height = 11
psfWidth3 = 68;
psfHeight3 = 11;

%% This function determines the psfTemplate based on psf dimensions
function [psfTemplate, box] = findPsfTemplate (frame, psfWidth, psfHeight, bubbleX, bubbleY)
    patchWidth = psfWidth*1.5;
    patchHeight = psfHeight*1.5;
    heightDiff = (patchHeight - psfHeight)/2;
    widthDiff = (patchWidth - psfWidth)/2;
    x = bubbleX - (psfWidth)/2 - widthDiff;
    y = bubbleY - (psfHeight)/2 - heightDiff;
    w = patchWidth;
    h = patchHeight;
    box = [x,y, w, h];
    % box = [4.685000000000000e+02, 1.115500000000000e+03, patchWidth, patchHeight]
    % box = [4.685000000000000e+02, 1.115500000000000e+03, patchWidth, patchHeight]
    % x = round(box(1));
    % y = round(box(2));
    % w = round(box(3));
    % h = round(box(4));
    y
    h
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

function [centroids] = locateCentroids (refFrame)
    refFrame = im2gray(refFrame);
    threshold = prctile(refFrame(:), 99);  % top 1% intensities, arbitrarily selected to include starburst
    regionalMax = imregionalmax(refFrame);
    bw = regionalMax & (refFrame > threshold);
    cc = bwconncomp(bw, 8);
    centroids = regionprops(cc, refFrame, 'WeightedCentroid');
end
% Locate centroids in the reference frame
% refFrame = refFrame1; % Use the line below to pick an appropriate frame
refFrame = read(simvid, 70);
centroids = locateCentroids(refFrame1);
figure;
imshow(refFrame)
impixelinfo;
title('Find Bubble Centres')








%lsqcurvefit



%% Discretation calculations, ignore for the most part
% ≥0.5λ in the lateral direction , ≥ 0.25λ in the axial direction.
fc = 2.0*1e6;
c = 1540;
lambda = c /fc; % λ = c / f
pixelSize = 3.3879*1e-5;
lambda/4 > pixelSize; % This returns true, so discretation is enough
