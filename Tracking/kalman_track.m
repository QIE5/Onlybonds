%% run_kalman_tracking.m
% Purpose: run "localise -> Kalman predict -> Hungarian assign -> correct" on simulation.mp4
% Dependencies:
%   - Computer Vision Toolbox: configureKalmanFilter, assignDetectionsToTracks
%   - Image Processing Toolbox: bwconncomp, regionprops, imregionalmax
% % % Files:
%   - localisationFunc.m must be in the same folder / on MATLAB path

clear; clc; close all;

%% --------------------- 0) Input video ---------------------
videoFile = "simulation.mp4";
assert(exist(videoFile,'file')==2, "Video file not found: %s", videoFile);

vid = VideoReader(videoFile);
fps = vid.FrameRate;

fprintf("Using video: %s | fps=%.3f | duration=%.2fs\n", videoFile, fps, vid.Duration);

maxFramesToProcess = 1500;   % You can reduce this for quick testing
doVisualise = true;

%% --------------------- 1) Localisation parameters (aligned with your existing code) ---------------------
loc.param = struct();
loc.param.threshold = 99.5;   % Same as in your localisation.m (percentile threshold)

loc.psf = struct();
loc.psf.height = 33;
loc.psf.width  = 183;

%% --------------------- 2) Tracking parameters (fps gating + Kalman noises) ---------------------
trk = struct();

% You used pixelSize (unit: m/px) in centroid.m; reuse it here
trk.pixelSize_m = 3.3879e-5;  % TODO: If you change the dataset later, recompute this

% "fps method": maxDist = vMax / fps / pixelSize
trk.vMax_mps = 0.10;          % TODO: Tune this according to flow / simulation settings
trk.maxDist_px = (trk.vMax_mps / fps) / trk.pixelSize_m;
trk.costOfNonAssignment = trk.maxDist_px;  % Cost is pixel distance, so the threshold is in pixels too

fprintf("max feasible distance ~= %.2f px/frame (vMax=%.3f m/s)\n", trk.maxDist_px, trk.vMax_mps);

% Track management
trk.invisibleForTooLong = 5;  % Delete tracks after this many consecutive missed frames
trk.minAge = 3;               % Very short tracks are usually unreliable
trk.minVisibility = 0.60;     % Delete tracks if visibility ratio is below this

% ✅ configureKalmanFilter expects these to be scalars (your previous error came from here)
trk.initEstimateError = 10;   % Initial estimate error (in pixel^2 sense)
trk.motionNoise       = 1;    % Motion noise (larger = allow stronger acceleration/maneuver)
trk.measurementNoise  = 4;    % Measurement noise (larger = trust detections less; smoother tracks)

%% --------------------- 3) Tracks container ---------------------
tracks = struct( ...
    'id', {}, ...
    'kalmanFilter', {}, ...
    'age', {}, ...
    'totalVisibleCount', {}, ...
    'consecutiveInvisibleCount', {}, ...
    'history', {} ...   % Nx3: [frameIdx, x, z]
);

nextId = 1;
frameIdx = 0;

if doVisualise
    figure('Name','ULM Kalman Tracking','Color','w');
end

%% --------------------- 4) Main loop: per-frame localisation + tracking ---------------------
while hasFrame(vid) && frameIdx < maxFramesToProcess
    frameIdx = frameIdx + 1;
    frame = readFrame(vid);

    % --- (1) localisation: output Nx2 [x,z] ---
    dets = localisationFunc(frame, loc.param, loc.psf);
    if isempty(dets), dets = zeros(0,2); end

    % --- (2) predict next positions for all tracks ---
    nTracks = numel(tracks);
    predicted = zeros(nTracks, 2);
    for i = 1:nTracks
        predicted(i,:) = predict(tracks(i).kalmanFilter);
    end

    % --- (3) Hungarian assignment: cost = distance(predicted, detections) ---
    if nTracks > 0 && ~isempty(dets)
        cost = zeros(nTracks, size(dets,1));
        for i = 1:nTracks
            dif = dets - predicted(i,:);
            cost(i,:) = sqrt(sum(dif.^2, 2))'; % distance (pixels)
        end

        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, trk.costOfNonAssignment);
    else
        assignments = zeros(0,2);
        unassignedTracks = (1:nTracks)';
        unassignedDetections = (1:size(dets,1))';
    end

    % --- (4) update assigned tracks: correct + append history ---
    for k = 1:size(assignments,1)
        t = assignments(k,1);
        d = assignments(k,2);

        corrected = correct(tracks(t).kalmanFilter, dets(d,:));

        tracks(t).age = tracks(t).age + 1;
        tracks(t).totalVisibleCount = tracks(t).totalVisibleCount + 1;
        tracks(t).consecutiveInvisibleCount = 0;

        tracks(t).history(end+1,:) = [frameIdx, corrected(1), corrected(2)];
    end

    % --- (5) update unassigned tracks: only increase invisible count; save predicted into history (for plotting) ---
    for k = 1:numel(unassignedTracks)
        t = unassignedTracks(k);

        tracks(t).age = tracks(t).age + 1;
        tracks(t).consecutiveInvisibleCount = tracks(t).consecutiveInvisibleCount + 1;

        if nTracks > 0
            tracks(t).history(end+1,:) = [frameIdx, predicted(t,1), predicted(t,2)];
        end
    end

    % --- (6) create new tracks for unassigned detections ---
    for k = 1:numel(unassignedDetections)
        d = unassignedDetections(k);

        kf = configureKalmanFilter( ...
            'ConstantVelocity', ...
            dets(d,:), ...
            trk.initEstimateError, ...
            trk.motionNoise, ...
            trk.measurementNoise);

        newTrack.id = nextId;
        newTrack.kalmanFilter = kf;
        newTrack.age = 1;
        newTrack.totalVisibleCount = 1;
        newTrack.consecutiveInvisibleCount = 0;
        newTrack.history = [frameIdx, dets(d,1), dets(d,2)];

        tracks(end+1) = newTrack; %#ok<AGROW>
        nextId = nextId + 1;
    end

    % --- (7) delete lost tracks ---
    if ~isempty(tracks)
        ages = [tracks.age];
        totalVis = [tracks.totalVisibleCount];
        visibility = totalVis ./ ages;

        lost = ([tracks.consecutiveInvisibleCount] >= trk.invisibleForTooLong) | ...
               ((ages < trk.minAge) & (visibility < trk.minVisibility));

        tracks = tracks(~lost);
    end

    % --- (8) Visualisation (raw frame + detections + tracks) ---
    if doVisualise
        imshow(frame); hold on;

        if ~isempty(dets)
            plot(dets(:,1), dets(:,2), 'r+', 'MarkerSize', 6, 'LineWidth', 1);
        end

        for i = 1:numel(tracks)
            h = tracks(i).history;
            if size(h,1) >= 2
                plot(h(:,2), h(:,3), '-', 'LineWidth', 1);
                plot(h(end,2), h(end,3), 'o', 'MarkerSize', 4);
                text(h(end,2)+3, h(end,3), sprintf('%d', tracks(i).id), ...
                    'Color','y','FontSize',8);
            end
        end

        title(sprintf('Frame %d | tracks=%d | dets=%d', frameIdx, numel(tracks), size(dets,1)));
        hold off;
        drawnow;
    end
end

fprintf("Finished. Total tracks remaining = %d\n", numel(tracks));

%% --------------------- 5) Output: density map (simple ULM density map) ---------------------
allPts = [];
for i = 1:numel(tracks)
    allPts = [allPts; tracks(i).history(:,2:3)]; %#ok<AGROW> % [x,z]
end

if ~isempty(allPts)
    xEdges = 0.5 : 1 : (vid.Width + 0.5);
    zEdges = 0.5 : 1 : (vid.Height + 0.5);

    % histcounts2 expects (X,Y). Here we use (z,x) so the output matrix matches image row/col
    density = histcounts2(allPts(:,2), allPts(:,1), zEdges, xEdges);

    figure('Name','ULM Density Map','Color','w');
    imagesc(density);
    axis image;
    colormap hot;
    title('ULM density map (counts per pixel)');
    xlabel('x (pixel)'); ylabel('z (pixel)');
end

%% Optional: save tracks
save("tracks_kalman_hungarian.mat", "tracks", "trk", "loc");
