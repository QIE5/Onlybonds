function [velocity, location, weight] = estimateTrackVelocityLocal (tracks, ...
    localisedBubbleCoords, offsets, frameRate);

    % Preallocate
    numTracks = numel(tracks);
    velocity = NaN(numTracks, 2);
    location = NaN(numTracks, 2);
    weight = zeros(numTracks, 1);
    
    for i = 1:numTracks
        n_pts = numel(tracks(i).frameIndices);
        i_start = max(1, floor(0.1*n_pts));
        i_end   = min(n_pts, ceil(0.9*n_pts));
        use_range = i_start:i_end;
        
        if numel(use_range) >= 3
            frames_use = finishedTracks(i).frameIndices(use_range);
            global_use = finishedTracks(i).globalIndices(use_range);
            
            % Extract positions
            pos_use = zeros(numel(frames_use), 2);
            for j = 1:numel(frames_use)
                fr = frames_use(j);
                localIdx = global_use(j) - offsets(fr);
                pos_use(j,:) = points{fr}(localIdx, :);
            end
            
            % Convert to time
            t = frames_use / frame_rate;
            
            % Robust linear fit
            px = robustfit(t, pos_use(:,1));
            py = robustfit(t, pos_use(:,2));
            
            vx = px(2);
            vy = py(2);
            
            % Residual-based weight
            res_x = pos_use(:,1) - (px(1) + px(2)*t);
            res_y = pos_use(:,2) - (py(1) + py(2)*t);
            fit_error = mean(res_x.^2 + res_y.^2);
            
            % Store results
            velocity(i,:) = [vx, vy];
            
            % Location at temporal midpoint
            t_mid = mean(t);
            location(i,:) = [px(1) + px(2)*t_mid, py(1) + py(2)*t_mid];
            
            weight(i) = numel(use_range) / (fit_error + eps);
        end
    end
    
    % Remove invalid tracks
    valid = ~isnan(velocity(:,1));
    velocity = velocity(valid, :);
    location = location(valid, :);
    weight = weight(valid);
    
    % Normalize weights (optional)
    weight = weight / median(weight);
    
end