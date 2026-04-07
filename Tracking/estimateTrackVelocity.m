function [velocity, location, weight] = estimateTrackVelocity(trackStructs, points, offsets, frameRate)
% ESTIMATETRACKVELOCITYLOCAL Compute velocity for each track using linear regression
%
% INPUT:
%   trackStructs : struct array from kalmanHungarianTracker
%   points       : cell array, one cell per frame, each cell is N x 2 [x y]
%   offsets      : array of global index offsets per frame
%   frameRate    : acquisition frame rate in Hz
%
% OUTPUT:
%   velocity : n_tracks x 2 array [vx vy] in units/second
%   location : n_tracks x 2 array [x y], track centroid
%   weight   : n_tracks x 1 array, track length (number of detections)

    nTracks = numel(trackStructs);
    velocity = zeros(nTracks, 2);
    location = zeros(nTracks, 2);
    weight = zeros(nTracks, 1);
    
    for i = 1:nTracks
        frameIndices = trackStructs(i).frameIndices;
        globalIndices = trackStructs(i).globalIndices;
        nPts = numel(frameIndices);
        
        % Extract positions
        pos = zeros(nPts, 2);
        for j = 1:nPts
            fr = frameIndices(j);
            localIdx = globalIndices(j) - offsets(fr);
            pos(j, :) = points{fr}(localIdx, :);
        end
        
        % Convert frames to time
        t = frameIndices / frameRate;
        
        % Linear fit: position = a + b*t, where b is velocity
        px = polyfit(t, pos(:, 1), 1);
        py = polyfit(t, pos(:, 2), 1);
        
        % Extract velocity (slope)
        velocity(i, :) = [px(1), py(1)];
        
        % Location at temporal midpoint
        tMid = mean(t);
        location(i, :) = [polyval(px, tMid), polyval(py, tMid)];
        
        % Weight by track length
        weight(i) = nPts;
    end
end