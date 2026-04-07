function [tracks, adjacency_tracks, A, trackStructs, offsets] = kalmanHungarianTracker1(points, varargin)
% KALMANHUNGARIANTRACKER
% Kalman prediction + Hungarian assignment tracker for ULM pipeline.
%
% INPUT:
%   points : cell array, one cell per frame
%            each cell is N x 2, coordinates [x y]
%
% OUTPUT:
%   tracks            : cell array, each track is n_frames x 1
%   adjacency_tracks  : cell array of global point indices
%   A                 : sparse adjacency matrix over concatenated points
%
% Optional key/value pairs:
%   'Debug'                  : false by default
%   'MaxGapClosing'          : 3 by default
%   'MaxLinkingDistance'     : 30 by default
%   'ProcessNoise'           : 10 by default
%   'MeasurementNoise'       : 4 by default
%   'InitialVelocityVariance': 100 by default
%   'MinTrackLength'         : 1 by default

    %% Parse arguments
    p = inputParser;
    p.addParameter('Debug', false, @(x) islogical(x) || isnumeric(x));
    p.addParameter('MaxGapClosing', 3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
    p.addParameter('MaxLinkingDistance', 30, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('ProcessNoise', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('MeasurementNoise', 4, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('InitialVelocityVariance', 100, @(x) isnumeric(x) && isscalar(x) && x > 0);
    p.addParameter('MinTrackLength', 1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    p.parse(varargin{:});

    debug                     = logical(p.Results.Debug);
    max_gap_closing           = p.Results.MaxGapClosing;
    max_linking_distance      = p.Results.MaxLinkingDistance;
    process_noise             = p.Results.ProcessNoise;
    measurement_noise         = p.Results.MeasurementNoise;
    initial_velocity_variance = p.Results.InitialVelocityVariance;
    min_track_length          = p.Results.MinTrackLength;

    %% Basic info
    n_slices = numel(points);

    if n_slices == 0
        tracks = cell(0,1);
        adjacency_tracks = cell(0,1);
        A = sparse(0,0);
        return;
    end

    n_cells = zeros(n_slices, 1);
    for i = 1:n_slices
        if isempty(points{i})
            n_cells(i) = 0;
        else
            n_cells(i) = size(points{i}, 1);
        end
    end

    n_total_cells = sum(n_cells);

    offsets = zeros(n_slices, 1);
    if n_slices > 1
        offsets(2:end) = cumsum(n_cells(1:end-1));
    end

    %% Kalman model
    % State = [x; y; vx; vy]
    % Here dt is treated as 1 frame.
    F = [1 0 1 0;
         0 1 0 1;
         0 0 1 0;
         0 0 0 1];

    H = [1 0 0 0;
         0 1 0 0];

    Q = process_noise * [0.25 0    0.5  0;
                         0    0.25 0    0.5;
                         0.5  0    1    0;
                         0    0.5  0    1];

    R = measurement_noise * eye(2);

    %% Track containers
    emptyTrack = struct( ...
        'id', [], ...
        'state', [], ...
        'P', [], ...
        'age', 0, ...
        'totalVisibleCount', 0, ...
        'consecutiveInvisibleCount', 0, ...
        'globalIndices', [], ...
        'frameIndices', [] ...
    );

    activeTracks   = repmat(emptyTrack, 0, 1);
    finishedTracks = repmat(emptyTrack, 0, 1);
    nextTrackID = 1;

    %% Main loop
    for frameIdx = 1:n_slices

        detections = points{frameIdx};

        if isempty(detections)
            detections = zeros(0,2);
        else
            detections = double(detections);
            if size(detections,2) > 2
                detections = detections(:,1:2);
            elseif size(detections,2) < 2
                error('Each points{t} must be an N x 2 array.');
            end
        end

        if debug
            fprintf('Frame %d/%d | detections = %d | active before predict = %d\n', ...
                frameIdx, n_slices, size(detections,1), numel(activeTracks));
        end

        %% Predict all active tracks
        nActive = numel(activeTracks);
        predictedPositions = zeros(nActive, 2);

        for i = 1:nActive
            [activeTracks(i).state, activeTracks(i).P] = ...
                predictKalman(activeTracks(i).state, activeTracks(i).P, F, Q);

            activeTracks(i).age = activeTracks(i).age + 1;
            predictedPositions(i,:) = activeTracks(i).state(1:2).';
        end

        %% Hungarian assignment between predicted positions and current detections
        [assignedPairs, unassignedTrackIdx, unassignedDetIdx] = ...
            hungarianAssign(predictedPositions, detections, max_linking_distance);

        %% Update assigned tracks
        for k = 1:size(assignedPairs,1)
            tIdx = assignedPairs(k,1);
            dIdx = assignedPairs(k,2);

            z = detections(dIdx,:).';

            [activeTracks(tIdx).state, activeTracks(tIdx).P] = ...
                correctKalman(activeTracks(tIdx).state, activeTracks(tIdx).P, z, H, R);

            activeTracks(tIdx).totalVisibleCount = activeTracks(tIdx).totalVisibleCount + 1;
            activeTracks(tIdx).consecutiveInvisibleCount = 0;

            globalIdx = offsets(frameIdx) + dIdx;
            activeTracks(tIdx).globalIndices(end+1,1) = globalIdx;
            activeTracks(tIdx).frameIndices(end+1,1) = frameIdx;
        end

        %% Mark unmatched tracks as invisible
        for k = 1:numel(unassignedTrackIdx)
            idx = unassignedTrackIdx(k);
            activeTracks(idx).consecutiveInvisibleCount = ...
                activeTracks(idx).consecutiveInvisibleCount + 1;
        end

        %% Remove dead tracks
        toRemove = false(numel(activeTracks),1);
        for i = 1:numel(activeTracks)
            if activeTracks(i).consecutiveInvisibleCount > max_gap_closing
                finishedTracks(end+1,1) = activeTracks(i); %#ok<AGROW>
                toRemove(i) = true;
            end
        end
        activeTracks(toRemove) = [];

        %% Create new tracks from unmatched detections
        for k = 1:numel(unassignedDetIdx)
            dIdx = unassignedDetIdx(k);

            x0 = [detections(dIdx,1); detections(dIdx,2); 0; 0];
            P0 = diag([measurement_noise, measurement_noise, ...
                       initial_velocity_variance, initial_velocity_variance]);

            tr = emptyTrack;
            tr.id = nextTrackID;
            tr.state = x0;
            tr.P = P0;
            tr.age = 1;
            tr.totalVisibleCount = 1;
            tr.consecutiveInvisibleCount = 0;
            tr.globalIndices = offsets(frameIdx) + dIdx;
            tr.frameIndices = frameIdx;

            activeTracks(end+1,1) = tr; %#ok<AGROW>
            nextTrackID = nextTrackID + 1;
        end

        if debug
            fprintf('Frame %d/%d | assigned = %d | active after update = %d | finished = %d\n', ...
                frameIdx, n_slices, size(assignedPairs,1), numel(activeTracks), numel(finishedTracks));
        end
    end

    %% Flush remaining active tracks
    if ~isempty(activeTracks)
        finishedTracks = [finishedTracks; activeTracks];
    end

    %% Filter short tracks
    if isempty(finishedTracks)
        tracks = cell(0,1);
        adjacency_tracks = cell(0,1);
        A = sparse(n_total_cells, n_total_cells);
        return;
    end

    keepMask = false(numel(finishedTracks),1);
    for i = 1:numel(finishedTracks)
        keepMask(i) = numel(finishedTracks(i).globalIndices) >= min_track_length;
    end
    finishedTracks = finishedTracks(keepMask);

    if isempty(finishedTracks)
        tracks = cell(0,1);
        adjacency_tracks = cell(0,1);
        A = sparse(n_total_cells, n_total_cells);
        return;
    end

    %% Sort tracks by first frame, then first global index
    sortKey = zeros(numel(finishedTracks), 2);
    for i = 1:numel(finishedTracks)
        sortKey(i,1) = finishedTracks(i).frameIndices(1);
        sortKey(i,2) = finishedTracks(i).globalIndices(1);
    end
    [~, order] = sortrows(sortKey, [1 2]);
    finishedTracks = finishedTracks(order);

    %% Build outputs
    n_tracks = numel(finishedTracks);
    tracks = cell(n_tracks,1);
    adjacency_tracks = cell(n_tracks,1);

    rows = [];
    cols = [];

    for i = 1:n_tracks
        adj = finishedTracks(i).globalIndices(:);
        frs = finishedTracks(i).frameIndices(:);

        adjacency_tracks{i} = adj;

        tr = NaN(n_slices,1);
        for j = 1:numel(adj)
            fr = frs(j);
            tr(fr) = adj(j) - offsets(fr);
        end
        tracks{i} = tr;

        if numel(adj) >= 2
            rows = [rows; adj(1:end-1)]; %#ok<AGROW>
            cols = [cols; adj(2:end)]; %#ok<AGROW>
        end
    end

    if isempty(rows)
        A = sparse(n_total_cells, n_total_cells);
    else
        A = sparse(rows, cols, 1, n_total_cells, n_total_cells);
    end
    %% Additional output variables for velocity calculations:
    trackStructs = finishedTracks;

end

%% ===== Kalman predict =====
function [xPred, PPred] = predictKalman(x, P, F, Q)
    xPred = F * x;
    PPred = F * P * F' + Q;
    PPred = (PPred + PPred') / 2;
end

%% ===== Kalman correct =====
function [xCorr, PCorr] = correctKalman(xPred, PPred, z, H, R)
    S = H * PPred * H' + R;
    K = (PPred * H') / S;

    innovation = z - H * xPred;
    xCorr = xPred + K * innovation;

    I = eye(size(PPred));
    PCorr = (I - K * H) * PPred * (I - K * H)' + K * R * K';
    PCorr = (PCorr + PCorr') / 2;
end

%% ===== Hungarian assignment wrapper =====
function [assignedPairs, unassignedTrackIdx, unassignedDetIdx] = ...
    hungarianAssign(predictedPositions, detections, maxDist)

    nT = size(predictedPositions,1);
    nD = size(detections,1);

    if nT == 0 && nD == 0
        assignedPairs = zeros(0,2);
        unassignedTrackIdx = zeros(0,1);
        unassignedDetIdx = zeros(0,1);
        return;
    elseif nT == 0
        assignedPairs = zeros(0,2);
        unassignedTrackIdx = zeros(0,1);
        unassignedDetIdx = (1:nD).';
        return;
    elseif nD == 0
        assignedPairs = zeros(0,2);
        unassignedTrackIdx = (1:nT).';
        unassignedDetIdx = zeros(0,1);
        return;
    end

    [target_indices, ~, unassigned_targets] = ...
        hungarianLinkerLocal(predictedPositions, detections, maxDist);

    matchedMask = target_indices ~= -1;
    matchedTrackIdx = find(matchedMask);
    matchedDetIdx = target_indices(matchedMask);

    assignedPairs = [matchedTrackIdx(:), matchedDetIdx(:)];
    unassignedTrackIdx = find(~matchedMask);
    unassignedDetIdx = unassigned_targets(:);
end

%% ===== Local Hungarian linker =====
function [target_indices, target_distances, unassigned_targets, total_cost] = ...
    hungarianLinkerLocal(source, target, max_distance)

    if nargin < 3
        max_distance = Inf;
    end

    n_source_points = size(source, 1);
    n_target_points = size(target, 1);

    %% Edge cases
    if n_source_points == 0
        target_indices = zeros(0,1);
        target_distances = zeros(0,1);
        unassigned_targets = (1:n_target_points).';
        total_cost = 0;
        return;
    end

    if n_target_points == 0
        target_indices = -1 * ones(n_source_points, 1);
        target_distances = NaN(n_source_points, 1);
        unassigned_targets = zeros(0,1);
        total_cost = Inf;
        return;
    end

    %% Build distance matrix
    D = NaN(n_source_points, n_target_points);

    for i = 1:n_source_points
        current_point = source(i, :);
        diff_coords = target - repmat(current_point, n_target_points, 1);
        square_dist = sum(diff_coords.^2, 2);
        D(i, :) = square_dist;
    end

    %% Apply max linking distance
    D(D > max_distance * max_distance) = Inf;

    %% If nothing is matchable, do not call munkres
    if all(isinf(D(:)))
        target_indices = -1 * ones(n_source_points, 1);
        target_distances = NaN(n_source_points, 1);
        unassigned_targets = (1:n_target_points).';
        total_cost = Inf;
        return;
    end

    %% Check munkres availability
    if exist('munkres', 'file') ~= 2
        error(['munkres.m not found. The Hungarian assignment in this tracker ' ...
               'requires the same munkres implementation used by the original simple tracker.']);
    end

    %% Hungarian solve
    [target_indices, total_cost] = munkres(D);
    target_indices(target_indices == 0) = -1;

    %% Distances
    target_distances = NaN(numel(target_indices), 1);
    for i = 1:numel(target_indices)
        if target_indices(i) < 0
            continue
        end
        target_distances(i) = sqrt(D(i, target_indices(i)));
    end

    %% Unassigned targets
    unassigned_targets = setdiff((1:n_target_points).', target_indices(target_indices > 0));
end