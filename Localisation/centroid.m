%% Discretation calculations
% ≥0.5λ in the lateral direction , ≥ 0.25λ in the axial direction.
fc = 2.0*1e6;
c = 1540;
lambda = c /fc; % λ = c / f
pixelSize = 3.3879*1e-5;
lambda/4 > pixelSize; % This returns true, so discretation is enough

%%
function [zc,xc] = localizeRadialSymmetry(I, fwhmz, fwhmx)
%% function [zc,xc] = localizeRadialSymmetry(I,fwhmz,fwhmx)
% Performs localization using radial symmetry properties
%
% Created by Baptiste Heiles on 05/09/18
% Inspired from Raghuveer Parthasarathy, The University of Oregon
%
% DATE 2020.07.22 - VERSION 1.1
% AUTHORS: Arthur Chavignon, Baptiste Heiles, Vincent Hingot. CNRS, Sorbonne Universite, INSERM.
% Laboratoire d'Imagerie Biomedicale, Team PPM. 15 rue de l'Ecole de Medecine, 75006, Paris
% Code Available under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (see https://creativecommons.org/licenses/by-nc-sa/4.0/)
% ACADEMIC REFERENCES TO BE CITED
% Details of the code in the article by Heiles, Chavignon, Hingot, Lopez, Teston and Couture.  
% Performance benchmarking of microbubble-localization algorithms for ultrasound localization microscopy, Nature Biomedical Engineering, 2021.
% General description of super-resolution in: Couture et al., Ultrasound localization microscopy and super-resolution: A state of the art, IEEE UFFC 2018
%
% Calculates the center of a 2D intensity distribution.
% Method: The gradient of a function that has perfect radial symmetry will
% point towards the origin. Thus we take the local gradient and construct
% lines through any point with orientation parallel to the local gradient.
% The origin is the point that will minimize the distance between itself
% and all such lines.
%
% INPUTS:
%       - I : 2D intensity distribution
%        Size need not be an odd number of pixels along each dimension
%       - fwhmz, fwhmx : full width at half maximum in direction z and x (unused)
% OUTPUTS:
%       - [zc, xc] : the center of radial symmetry,
%            px, from center

%% Number of grid points
[Nz,Nx] = size(I);

%% Radial symmetry algorithm

% grid coordinates are -n:n, where Nz (or Nx) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5. Note that z increases "downward"
zm_onerow = (-(Nz-1)/2.0+0.5:(Nz-1)/2.0-0.5)';
zm = zm_onerow(:,ones(Nx-1, 1));
xm_onecol = (-(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5);
xm = xm_onecol(ones(Nz-1, 1),:);

% Calculate derivatives along 45-degree shifted coordinates (u and v) Please refer to Appendix 2 of the publication attached to this code for basis definition
dIdu = I(1:Nz-1,2:Nx)-I(2:Nz,1:Nx-1);% Gradient along the u vector
dIdv = I(1:Nz-1,1:Nx-1)-I(2:Nz,2:Nx);% Gradient along the v vector

% Smoothing the gradient of the I window
h = ones(3)/9;
fdu = conv2(dIdu, h, 'same');% Convolution of the gradient with a simple averaging filter
fdv = conv2(dIdv, h, 'same');
dImag2 = fdu.*fdu + fdv.*fdv; % Squared gradient magnitude

% Slope of the gradient . Please refer to appendix 2 of the publication attached to this code for basis/orientation
m = -(fdv + fdu) ./ (fdu-fdv);

% Check if m is NaN (which can happen when fdu=fdv). In this case, replace with the un-smoothed gradient.
NNanm = sum(isnan(m(:)));
if NNanm > 0
    unsmoothm = (dIdv + dIdu) ./ (dIdu-dIdv);
    m(isnan(m))=unsmoothm(isnan(m));
end

% If it's still NaN, replace with zero and we'll deal with this later
NNanm = sum(isnan(m(:)));
if NNanm > 0
    m(isnan(m))=0;
end

% Check if m is inf (which can happen when fdu=fdv).
try
    m(isinf(m))=10*max(m(~isinf(m)));
catch
    % Replace m with the unsmoothed gradient
    m = (dIdv + dIdu) ./ (dIdu-dIdv);
end

% Calculate the z intercept of the line of slope m that goes through each grid midpoint
b = zm - m.*xm;

% Weight the intensity by square of gradient magnitude and inverse
% distance to gradient intensity centroid. This will increase the intensity of areas close to the initial guess
sdI2 = sum(dImag2(:));
zcentroid = sum(sum(dImag2.*zm))/sdI2;% Initial guess of the centroid in z
xcentroid = sum(sum(dImag2.*xm))/sdI2;% Initial guess of the centroid in x
w  = dImag2./sqrt((zm-zcentroid).*(zm-zcentroid)+(xm-xcentroid).*(xm-xcentroid));

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc) such that lines y = mx+b have
% the minimal total distance^2 to the origin:
% See function lsradialcenterfit (below)
[zc,xc] = lsradialcenterfit(m, b, w);

end

% We'll code the least square solution function separately as we could find the solution with another implementation
function [zc,xc] = lsradialcenterfit(m, b, w)
    % least squares solution to determine the radial symmetry center

    % inputs m, b, w are defined on a grid
    % w are the weights for each point
    wm2p1 = w./(m.*m+1);
    sw  = sum(sum(wm2p1));
    smmw = sum(sum(m.*m.*wm2p1));
    smw  = sum(sum(m.*wm2p1));
    smbw = sum(sum(m.*b.*wm2p1));
    sbw  = sum(sum(b.*wm2p1));
    det = smw*smw - smmw*sw;
    xc = (smbw*sw - smw*sbw)/det;    % relative to image center
    zc = (smbw*smw - smmw*sbw)/det; % relative to image center

end



%% Code for producing video with localied MBs
simvid = VideoReader('simulation.mp4');
% localisedvid = VideoWriter("localised");
% localisedvid.FrameRate = simvid.FrameRate;
% open(localisedvid);
% figure('Visible','off');

numFrames = simvid.NumFrames;
centroidCell = cell(numFrames, 1);
frame1 = readFrame(simvid);
figure;
imshow(frame1);
hold on
plot(417.4997 ,1.2695e+03,'b*');
% The lines below are drawn through the centre of the selected MB.
% They represent the dimensions of the MB
plot ([417.4997 - 90 ,417.4997 + 90], [ 1.2695e+03,1.2695e+03], 'r-'); % PSF width = 180
plot ([417.4997,417.4997], [ 1.2695e+03 + 15,1.2695e+03 - 15], 'g-'); % heaight = 30
hold off
patchSize = 180 + 3;
%% Trying radial
% Parameters
fwhmz = 180;  % arbitrary, adjust to your bubble size
fwhmx = 30;

% FWHM vectors (range around each peak)
vectfwhmz = -round(fwhmz/2):round(fwhmz/2);  % e.g., [-1, 0, 1]
vectfwhmx = -round(fwhmx/2):round(fwhmx/2);

% Margins for bounds checking
marginZ = round(fwhmz/2);
marginX = round(fwhmx/2);
%%
numFrames = 2;
for n = 1: numFrames;
    frame = read(simvid, n);
    frame = im2gray(frame);
    frame = double(frame);
    threshold = prctile(frame(:), 99.5);  % top 0.5% intensities, arbitrarily selected to include starburst
    regionalMax = imregionalmax(frame);
    % bw = frame > threshold;
    bw = regionalMax & (frame > threshold);
    % [xPeaks, yPeaks] = find(bw)
    [zPeaks, xPeaks] = find(bw);
    
    cc = bwconncomp(bw, 8);
    nPeaks = cc.NumObjects;
finalZ = zeros(nPeaks, 1);
finalX = zeros(nPeaks, 1);

for iPeak = 1:nPeaks
    pixelIdx = cc.PixelIdxList{iPeak};
    [zAll, xAll] = ind2sub(size(frame), pixelIdx);
    z0 = round(mean(zAll));  % simple mean for initial ROI center
    x0 = round(mean(xAll));
    IntensityRoi = frame(z0 + vectfwhmz, x0 + vectfwhmx);
    [Zc, Xc] = localizeRadialSymmetry(IntensityRoi, fwhmz, fwhmx);
    finalZ(iPeak) = Zc + z0;
    finalX(iPeak) = Xc + x0;
end
    stats = regionprops(cc, frame, 'WeightedCentroid');
    % 
    cat(1,stats.WeightedCentroid)
    % centroidCell{n} = cat(1,stats.WeightedCentroid)
    centroidCell{n} = cat(1,[finalX,finalZ]);
    imshow(frame); hold on
    plot(centroidCell{n}(:,1), centroidCell{n}(:,2), 'r+', 'MarkerSize',6)
    hold off
    % 
    % F = getframe(gca);
    % writeVideo(localisedvid, F);
end
%%
centroids = centroidCell{1}

figure;
imshow(frame1);
hold on

plot(centroids(:,1),centroids(:,2),'b*');
hold off
title('Frame 1');

% close(localisedvid);
% close(gcf);


