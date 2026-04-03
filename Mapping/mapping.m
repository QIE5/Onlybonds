function [SR_img, BW] = mapping(frame, localisedBubbleCoords, tracks, adjacency_tracks, A)
    % === 1. CONCATENATE RAW LOCALIZATIONS ===
    allPointsRaw = vertcat(localisedBubbleCoords{:});   % [x z]
    
    % === 2. BUILD INTERPOLATED POINTS FROM TRACKS ===
    allPoints = [];
    
    for k = 1:numel(adjacency_tracks)
        idx = adjacency_tracks{k};
        
        if numel(idx) < 2, continue; end
        
        P = allPointsRaw(idx, :);
        
        s  = 1:size(P,1);
        si = linspace(1, size(P,1), 10*numel(s));
        
        xi = interp1(s, P(:,1), si);
        zi = interp1(s, P(:,2), si);
        
        allPoints = [allPoints; xi(:), zi(:)];
    end
    
    % === 3. MAPPING ===
    [H,W] = size(frame);
    
    xEdges = linspace(min(allPoints(:,1)), max(allPoints(:,1)), W+1);
    zEdges = linspace(min(allPoints(:,2)), max(allPoints(:,2)), H+1);
    
    SR_count = histcounts2(allPoints(:,2), allPoints(:,1), zEdges, xEdges);
    
    % === 4. PSF RENDERING ===
    SR_img = imgaussfilt(SR_count, 1.5);
    
    % === 5. VESSEL SEGMENTATION ===
    nz  = SR_img(SR_img > 0);
    thr = prctile(nz, 20);
    BW  = SR_img > thr;
    
    BW = imclose(BW, strel('disk',5));
    BW = imfill(BW, 'holes');
    BW = bwareaopen(BW, 50);
end


% === 6. DISPLAY ===
% figure;
% subplot(1,2,1);
% imagesc(log(SR_img + 1));
% axis image off; colormap hot;
% title('ULM Density Map');
% 
% subplot(1,2,2);
% imshow(BW);
% title('Reconstructed Vessel');

% Can overlay on B-mode image
