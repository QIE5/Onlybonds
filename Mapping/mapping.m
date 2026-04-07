function [srImg, velImg] = mapping(frame, localisedBubbleCoords, adjacency_tracks, velocity, location, weight, frameRate)
% MAPPING Generate super-resolved structure and velocity images
%
% INPUT:
%   frame                 : original B-mode frame (for dimensions)
%   localisedBubbleCoords : cell array of localizations [x z]
%   adjacency_tracks      : cell array of global indices per track
%   velocity              : Nx2 array [vx, vz]
%   location              : Nx2 array [x, z] positions
%   weight                : Nx1 array of weights
%   frameRate             : acquisition frame rate in Hz
%
% OUTPUT:
%   srImg  : super-resolved density image
%   velImg : struct with .vx, .vy, .speed, .direction, .validMask

    scale = 1;
    filterSigma = 1.5;  % Gaussian smoothing sigma in pixels

    [H, W] = size(frame);
    xEdges = linspace(1, W, scale*W + 1);
    zEdges = linspace(1, H, scale*H + 1);
    
    rows = numel(zEdges) - 1;
    cols = numel(xEdges) - 1;

    % ================================================================
    % STRUCTURE IMAGE
    % ================================================================

    allPointsRaw = vertcat(localisedBubbleCoords{:});

    totalPts = 0;
    for k = 1:numel(adjacency_tracks)
        idx = adjacency_tracks{k};
        if numel(idx) < 2, continue; end
        totalPts = totalPts + 10 * numel(idx);
    end

    allPoints = zeros(totalPts, 2);
    ptr = 1;
    for k = 1:numel(adjacency_tracks)
        idx = adjacency_tracks{k};
        if numel(idx) < 2, continue; end
        P = allPointsRaw(idx, :);
        s  = 1:size(P, 1);
        si = linspace(1, size(P, 1), 10*numel(s));
        xi = interp1(s, P(:,1), si);
        zi = interp1(s, P(:,2), si);
        n = numel(xi);
        allPoints(ptr:ptr+n-1, :) = [xi(:), zi(:)];
        ptr = ptr + n;
    end
    allPoints = allPoints(1:ptr-1, :);

    srImg = histcounts2(allPoints(:,2), allPoints(:,1), zEdges, xEdges);
    srImg = imgaussfilt(srImg, filterSigma);

    % ================================================================
    % VELOCITY IMAGE
    % ================================================================
    if nargout > 1
        weight = weight(:);
        n = size(location, 1);
        radius = 30;

        % ---- Step 1: Velocity smoothing (Eqs 3-5) via KD-tree ----
        neighborIdx = rangesearch(location, location, radius);

        vxSmoothed = zeros(n, 1);
        vySmoothed = zeros(n, 1);
        r2 = radius^2;

        for i = 1:n
            idx = neighborIdx{i};
            if isempty(idx)
                vxSmoothed(i) = velocity(i,1);
                vySmoothed(i) = velocity(i,2);
                continue;
            end

            dx = location(idx,1) - location(i,1);
            dz = location(idx,2) - location(i,2);
            d2 = dx.^2 + dz.^2;

            wi = exp(-d2 / r2) .* weight(idx);
            Z = sum(wi);

            if Z > 0
                vxSmoothed(i) = sum(wi .* velocity(idx,1)) / Z;
                vySmoothed(i) = sum(wi .* velocity(idx,2)) / Z;
            else
                vxSmoothed(i) = velocity(i,1);
                vySmoothed(i) = velocity(i,2);
            end
        end

        % ---- Step 2: Bin onto same grid as srImg ----
        colIdx = discretize(location(:,1), xEdges);
        rowIdx = discretize(location(:,2), zEdges);

        valid = ~isnan(colIdx) & ~isnan(rowIdx);
        rowIdx = rowIdx(valid);
        colIdx = colIdx(valid);
        vxS = vxSmoothed(valid);
        vyS = vySmoothed(valid);
        wV  = weight(valid);

        vxAccum = accumarray([rowIdx, colIdx], wV .* vxS, [rows, cols], @sum, 0);
        vyAccum = accumarray([rowIdx, colIdx], wV .* vyS, [rows, cols], @sum, 0);
        wAccum  = accumarray([rowIdx, colIdx], wV,        [rows, cols], @sum, 0);

        % ---- Step 3: Same smoothing as srImg ----
        vxSmooth = imgaussfilt(vxAccum, filterSigma);
        vySmooth = imgaussfilt(vyAccum, filterSigma);
        wSmooth  = imgaussfilt(wAccum,  filterSigma);

        % ---- Step 4: Normalize where weight exists ----
        hasData = wSmooth > max(wSmooth(:)) * 1e-4;

        vxImage = zeros(rows, cols);
        vyImage = zeros(rows, cols);
        vxImage(hasData) = vxSmooth(hasData) ./ wSmooth(hasData);
        vyImage(hasData) = vySmooth(hasData) ./ wSmooth(hasData);

        % ---- Step 5: Use srImg as vessel mask ----
        vesselMask = srImg > max(srImg(:)) * 1e-3;

        % ---- Step 6: Fill gaps within vessel mask ----
        % Where vessels exist but velocity is missing, interpolate
        needFill = vesselMask & ~hasData;

        if any(needFill(:))
            vxImage = fillGaps(vxImage, hasData, needFill);
            vyImage = fillGaps(vyImage, hasData, needFill);
        end

        % Final mask: velocity exists everywhere vessels exist
        validMask = vesselMask;

        velImg = struct();
        velImg.vx = vxImage;
        velImg.vy = vyImage;
        velImg.speed = hypot(vxImage, vyImage);
        velImg.direction = atan2(vyImage, vxImage);
        velImg.count = wSmooth;
        velImg.validMask = validMask;
    else
        velImg = struct();
    end
end

function img = fillGaps(img, hasData, needFill)
% FILLGAPS Interpolate missing pixels from nearby data pixels
% Uses iterative nearest-neighbor diffusion within the gap region

    filled = hasData;
    maxIter = 50;
    kernel = [0 1 0; 1 0 1; 0 1 0];  % 4-connected neighbors

    for iter = 1:maxIter
        remaining = needFill & ~filled;
        if ~any(remaining(:))
            break;
        end

        % Sum of filled neighbor values
        sumNeighbors = conv2(img .* double(filled), kernel, 'same');
        % Count of filled neighbors
        countNeighbors = conv2(double(filled), kernel, 'same');

        % Pixels that can be filled this iteration
        canFill = remaining & countNeighbors > 0;

        if ~any(canFill(:))
            break;
        end

        img(canFill) = sumNeighbors(canFill) ./ countNeighbors(canFill);
        filled(canFill) = true;
    end
end