%% Discretation calculations
% ≥0.5λ in the lateral direction , ≥ 0.25λ in the axial direction.
fc = 2.0*1e6;
c = 1540;
lambda = c /fc; % λ = c / f
pixelSize = 3.3879*1e-5;
lambda/4 > pixelSize; % This returns true, so discretation is enough


% simvid = VideoReader('simulation.mp4');
simvid = VideoReader('simulation.mp4');

%% Locate 3 different microbubbles, 1 for each region in the frame. Run 
% the code using an appropriate frame. Then, use 
% the centre coordinates (from centroids) and approximate width and height 
% of the bubbles. Continue to re-run code until all 3 bubbles are 
% approximated. These properties are used to form a psfTemplate for each 
% region. Use the psfTemplates to run crossCorrelation.m



frame1 = read(simvid, 1);
frame30 = read(simvid, 30);
figure;
imshow(frame1);
impixelinfo;
hold on

%% Simulation.mp4 microbubbles
bubble1 = [508,1068]; % Frame 1
plot(bubble1(1), bubble1(2), 'b*');
bubble2 = [435,1489]; % Frame 30
% plot(bubble2(1), bubble2(2), 'b*');
bubble3 = [534, 2232]; % Frame 1
plot(bubble3(1), bubble3(2), 'b*');

% The lines below are drawn through the centre of the selected MB.
% They represent the dimensions of the MB
% Bubble 1
plot ([bubble1(1) - 100 ,bubble1(1) + 100], [bubble1(2), bubble1(2)], 'r-'); % PSF width = 200
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + 15,bubble1(2) - 15], 'g-'); % height = 30
psfWidth1 = 200;
psfHeight1 = 30;
% Bubble 2, frame 30
plot ([bubble2(1) - 85 ,bubble2(1) + 85], [bubble2(2), bubble2(2)], 'r-'); % PSF width = 170
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 14,bubble2(2) - 14], 'g-'); % height = 28
psfWidth2 = 170;
psfHeight2 = 28;
% Bubble 3
plot ([bubble3(1) - 70 ,bubble3(1) + 40], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 110
% plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 12,bubble2(2) - 12], 'g-'); % height = 24
psfWidth3 = 110;
psfHeight3 = 24;

%% Find the template based on psf dimensions
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
[psfTemplate1, box1] = findPsfTemplate (frame1, psfWidth1, psfHeight1, bubble1(1), bubble1(2));
[psfTemplate2, box2] = findPsfTemplate (frame30, psfWidth2, psfHeight2, bubble2(1), bubble2(2));
[psfTemplate3, box3] = findPsfTemplate (frame1, psfWidth3, psfHeight3, bubble3(1), bubble3(2));



rectangle('Position', box1, 'EdgeColor', 'g', 'LineWidth', 1);
rectangle('Position', box2, 'EdgeColor', 'g', 'LineWidth', 1);
rectangle('Position', box3, 'EdgeColor', 'g', 'LineWidth', 1);
hold off

%%
    frame1 = im2gray(frame1);
    threshold = prctile(frame1(:), 99.5);  % top 0.5% intensities, arbitrarily selected to include starburst
    regionalMax = imregionalmax(frame1);
    % bw = frame > threshold;
    bw = regionalMax & (frame1 > threshold);
    % [xPeaks, yPeaks] = find(bw)
    
    cc = bwconncomp(bw, 8);
    centroids = regionprops(cc, frame1, 'WeightedCentroid');



%%
% centroids = centroidCell{1};
% figure;
% imshow(frame1);
% hold on
% plot(centroids(:,1),centroids(:,2),'b*');
% hold off
% title('Frame 1');

% close(localisedvid);
% close(gcf);



%lsqcurvefit
