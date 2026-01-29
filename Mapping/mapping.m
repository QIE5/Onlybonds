allPointsRaw = vertcat(localisedBubbleCoords{:});   

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
[H,W] = size(frame);
xEdges = linspace(min(allPoints(:,1)), max(allPoints(:,1)), W+1);
zEdges = linspace(min(allPoints(:,2)), max(allPoints(:,2)), H+1);
SR_count = histcounts2(allPoints(:,2), allPoints(:,1), zEdges, xEdges);
figure;
subplot(1,2,1);
imagesc(log(SR_img + 1));
axis image off; colormap hot;
title('ULM Density Map');
subplot(1,2,2);
imshow(SR_img);
title('Reconstructed Vessel');

% vector direction, velocity