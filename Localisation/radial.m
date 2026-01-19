function MatTracking = ULM_localization2D_mesh(MatIn,ULM)

%% Init struct
if ~isfield(ULM,'LocMethod'),ULM.LocMethod = 'Radial';end
if ~isfield(ULM,'parameters'),ULM.parameters = struct();end
if ~isfield(ULM,'threshold'),ULM.threshold = 99.5;end  % percentile threshold

if strcmp(ULM.LocMethod,'Interp')
    if ~isfield(ULM.parameters,'InterpMethod')
        ULM.parameters.InterpMethod = 'spline';
    end
end

%% Fill parameters
LocMethod = ULM.LocMethod;
fwhmz = ULM.fwhm(2); 
fwhmx = ULM.fwhm(1);

%% Initializing variables
height = size(MatIn, 1);
width = size(MatIn, 2);
numberOfFrames = size(MatIn, 3);

MatIn = abs(double(MatIn));

MatInReduced = zeros(height, width, numberOfFrames);
redz = 3:height-2;
redx = 3:width-2;
MatInReduced(redz, redx, :) = MatIn(redz, redx, :);

%% Multi-bubble detection (REPLACES the old single-max code)
threshold = prctile(MatInReduced(:), ULM.threshold);
BW = imregionalmax(MatInReduced) & (MatInReduced > threshold);
[index_mask_z, index_mask_x, index_mask_frame] = ind2sub(size(MatInReduced), find(BW));

nBubbles = length(index_mask_z);
fprintf('Detected %d candidate bubbles\n', nBubbles);

%% Creating FWHM vectors
vectfwhmz = -round(fwhmz/2):round(fwhmz/2);
vectfwhmx = -round(fwhmx/2):round(fwhmx/2);

%% Margin check
marginZ = round(fwhmz/2);
marginX = round(fwhmx/2);

%% Localization
averageZc = nan(nBubbles, 1);
averageXc = nan(nBubbles, 1);
convergent_vals = nan(nBubbles, 1);  % track which localizations succeeded

for iscat = 1:nBubbles
    
    z0 = index_mask_z(iscat);
    x0 = index_mask_x(iscat);
    f0 = index_mask_frame(iscat);
    
    % Check bounds
    if z0 - marginZ < 1 || z0 + marginZ > height || ...
       x0 - marginX < 1 || x0 + marginX > width
        continue
    end
    
    IntensityRoi = MatIn(z0 + vectfwhmz, x0 + vectfwhmx, f0);
    
    switch LocMethod
        case 'Radial'
            [Zc, Xc, ~] = LocRadial(IntensityRoi, fwhmz, fwhmx);
        case 'WA'
            [Zc, Xc, ~] = LocWeightedAverage(IntensityRoi, vectfwhmz, vectfwhmx);
        case 'Interp'
            [Zc, Xc, ~] = LocInterp(IntensityRoi, ULM.parameters.InterpMethod, vectfwhmz, vectfwhmx);
        case 'CurveFitting'
            [Zc, Xc, ~] = curveFitting(IntensityRoi, vectfwhmz, vectfwhmx);
        case 'NoLocalization'
            [Zc, Xc, ~] = NoLocalization(IntensityRoi);
        otherwise
            error('Wrong LocMethod selected')
    end
    
    averageZc(iscat) = Zc + z0;
    averageXc(iscat) = Xc + x0;
    convergent_vals(iscat) = f0;
end

%% Remove failed localizations
validIdx = ~isnan(averageZc);
averageZc = averageZc(validIdx);
averageXc = averageXc(validIdx);
convergent_vals = convergent_vals(validIdx);

%% Output: [Z, X, Frame]
MatTracking = [averageZc, averageXc, convergent_vals];

fprintf('Successfully localized %d bubbles\n', size(MatTracking, 1));

end
%% Additional Functions

function sigma = ComputeSigmaScat(Iin,Zc,Xc)
[Nx,Nz] = size(Iin);
Isub = Iin - mean(Iin(:));
[px,pz] = meshgrid(1:Nx,1:Nz);
zoffset = pz - Zc+(Nz)/2.0;%BH xoffset = px - xc;
xoffset = px - Xc+(Nx)/2.0;%BH yoffset = py - yc;
r2 = zoffset.*zoffset + xoffset.*xoffset;
sigma = sqrt(sum(sum(Isub.*r2))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
end

function [zc,xc] = localizeRadialSymmetry(I,fwhmz,fwhmx)
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
function [Zc,Xc,sigma] = LocRadial(Iin,fwhm_z,fwhm_x)
%% function [Zc,Xc,sigma] = LocRadial(Iin,fwhm_z,fwhm_x)
[Zc,Xc] = localizeRadialSymmetry(Iin,fwhm_z,fwhm_x);
sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = NoLocalization(Iin)
%% function [Zc,Xc,sigma] = NoLocalization(Iin)
Xc = 0;
Zc = 0;

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
%%function [Zc,Xc,sigma] = LocWeightedAverage(Iin,vectfwhm_z,vectfwhm_x)
Zc = sum(sum(Iin.*vectfwhm_z',1),2)./sum(Iin(:));
Xc = sum(sum(Iin.*vectfwhm_x,1),2)./sum(Iin(:));

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)
%%function [Zc,Xc,sigma] = LocInterp(Iin,InterpMode,vectfwhm_z,vectfwhm_x)

Nz=size(Iin,1);Nx=size(Iin,2);
if strcmp(InterpMode,'spline')
    [X,Z] = meshgrid(1:Nx,1:Nz);
    [Xq,Zq] = meshgrid(linspace(1,Nx,Nx*10),linspace(1,Nz,Nz*10));%to avoid uneven shift
    In_interp = interp2(X,Z,Iin,Xq,Zq,InterpMode);
else
    %     if strcmp(InterpMode,'lanczos3')
    %         In_interp = imresize(Iin,10,{@mylanczos2,7});
    %     else
    In_interp = imresize(Iin,10,InterpMode);
    %     end
end
[~,tt] = max(In_interp(:));
% [tt]=find(abs(In_interp)==max(abs(In_interp(:))));
% if size(tt,1)>1
%     tt=tt(1,1);%arbitrary solve this
% end
[iz,ix,~]=ind2sub(size(In_interp),tt);

Zc = vectfwhm_z(1)-1 + 1 + iz./10 -.5 +.05;
Xc = vectfwhm_x(1)-1 + 1 + ix./10 -.5 +.05;

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
%% function [Zc,Xc,sigma] = curveFitting(Iin,vectfwhm_z,vectfwhm_x)
[meshX,meshZ] = meshgrid(vectfwhm_x,vectfwhm_z);
meshIn = cat(3,meshX,meshZ);

sigGauss_z = vectfwhm_z(end)*0+1;
sigGauss_x = vectfwhm_x(end)*0+1;

myGaussFunc = @(x_pos,mesh_pos)( exp(-(mesh_pos(:,:,1)-x_pos(1)).^2./(2*sigGauss_z^2) - (mesh_pos(:,:,2)-x_pos(2)).^2./(2*sigGauss_x^2)));
OPTIONS = optimoptions('lsqcurvefit','StepTolerance',.01,'MaxIterations',5,'Display','off');

% Gaussian Fitting
x_out = lsqcurvefit(myGaussFunc,[0 0],meshIn,double(Iin./max(Iin(:))),[],[],OPTIONS);

Zc = x_out(2);
Xc = x_out(1);

sigma = ComputeSigmaScat(Iin,Zc,Xc);
end

function f = mylanczos2(x)
f = (sin(pi*x) .* sin(pi*x/2) + eps) ./ ((pi^2 * x.^2 / 2) + eps);
f = f .* (abs(x) < 2);
end

simvid = VideoReader('simulation.mp4');
frame = readFrame(simvid);
frame1 = readFrame(simvid);
frame1 = im2gray(frame1);
frame1 = double(frame1);

% Create 3D array
MatIn = cat(3, frame1, frame1);   % 2362 × 894 × 2

% Verify
disp(size(MatIn))
disp(ndims(MatIn))

% Set parameters
ULM = struct();
ULM.LocMethod = 'CurveFitting';
ULM.fwhm = [3, 3];
ULM.threshold = 99.5;   % adjust if too many/few detections

% Run localization - use MatIn, not frame1!
MatTracking = ULM_localization2D_mesh(MatIn, ULM);

% Visualize frame 1 only
figure
imagesc(frame1)
colormap gray
axis image
hold on

% Plot only detections from frame 1
idx = MatTracking(:,3) == 1;
plot(MatTracking(idx,2), MatTracking(idx,1), 'r+', 'MarkerSize', 10, 'LineWidth', 2)
title(sprintf('Frame 1: %d bubbles detected', sum(idx)))


