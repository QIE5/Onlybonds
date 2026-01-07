
%LOCALISATION


simvid = VideoReader('simulation.mp4');


                                   
%bwhd = hdomes > (max(hdomes(:)) * binthresh); %binarize at threshold

function hdomes = morphReconstruct(vid, nfr, h)
% morphReconstruct - Performs morphological reconstruction on a video frame
% DOES NOT BINARIZE
% Syntax: [hdomes, cehd] = morphReconstruct(vid, h, nfr)
%
% Inputs:
%   vid  - Video object (VideoReader or similar)
%   h    - Morphological offset value. Used to obtain marker image.
%          Use 0.05 nominally. 
%          Lower h = higher sensitivity in peak detection
%   nfr  - Frame number to process
%
% Outputs:
%   hdomes - Morphologically reconstructed image (mask - dilation)
%   cehd   - Contrast-enhanced version of hdomes


frame = rgb2gray(read(vid, nfr));       %load in frame and convert to grayscale

mask = frame;                           %Create mask(original)                           
marker = mask * (1 - h);                %Create marker(shrunken OG)
dilation = imreconstruct(marker, mask); %Dilation. not sure how this works tbh

hdomes = mask - dilation; %subtract dilation from OG to get final img

end

function output = cehd(a)%contrast enhance hdomes. for ease
    output = a * (255 / max(a(:))); 
end





b = a > 0;

%imshow(cehd(a))
%imshow(b)

vid = simvid;
h = 0.05;
nfr = 0;


a = morphReconstruct(vid, 56, h);
imshow(cehd(a));

b = a > 0;



binthresh = 0.1 %
volth = 30; %volume threshold, this threshold feels arbitrary? 


%{
nfr = 500;
while nfr < 1700
    nfr = nfr+1
    
    imshow(bwf);
    a = morphReconstruct(vid, nfr, h); 
    imshow(rgb2gray(read(vid, nfr))); %display
    bwim = a > 0.1 * max(a(:)); %don't even know how this binarization works
    bwf = bwareaopen(bwim, volth); %removes blobs smaller than a certain area

    
    
end

%}


bwf = bwareaopen(b, volth); %removes blobs smaller than a certain area
imshow(bwf)