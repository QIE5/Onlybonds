function [correctedFrames] = motionCorrection(frames)
    % outputFile = 'stabilized_final.mp4';
    % vOut = VideoWriter(outputFile, 'MPEG-4');
    % vOut.FrameRate = v.FrameRate;
    % open(vOut);
    % 1. Setup Reference
    % refFrame = im2gray(readFrame(v));
    refFrame = frames(:,:, 1);
    [rows, cols] = size(refFrame);
    refFrameD = double(refFrame);
    
    % Create a Hanning window to reduce FFT edge artifacts
    win = hanning(rows) * hanning(cols)';
    
    % Pre-calculate Reference FFT
    fftRef = fft2(refFrameD .* win);
    
    % Store shifts for smoothing
    % totalFrames = floor(v.Duration * v.FrameRate);
    totalFrames = size(frames, 3);
    shifts = zeros(totalFrames, 2);
    
    % 2. Calculate Shifts (First Pass)
    % v.CurrentTime = 0;
    fprintf('Analyzing motion...\n');
    for i = 1:totalFrames
        % if ~hasFrame(v), break; end
        % curr = double(im2gray(readFrame(v)));
        curr = double(im2gray(frames(:,:,i)));
        
        % Phase Correlation
        fftCurr = fft2(curr .* win);
        % Cross-power spectrum
        R = (fftRef .* conj(fftCurr)) ./ (abs(fftRef .* conj(fftCurr)) + eps);
        shiftMap = real(ifft2(R));
        
        [maxVal, maxIdx] = max(shiftMap(:));
        [rp, cp] = ind2sub([rows, cols], maxIdx);
        
        % --- SUB-PIXEL ESTIMATION (Parabolic Fit) ---
        % Check neighbors to interpolate the true peak location
        rRange = mod(rp-2:rp, rows) + 1;
        cRange = mod(cp-2:cp, cols) + 1;
        patch = shiftMap(rRange, cRange);
        
        % Parabolic interpolation for Sub-pixel accuracy
        % (Prevents the "jitter" caused by rounding to the nearest pixel)
        y = patch(:, 2);
        x = patch(2, :);
        
        % Sub-pixel offset in rows
        if y(1) ~= y(3) % Avoid division by zero
            dr_sub = (y(3) - y(1)) / (2 * (2 * y(2) - y(1) - y(3)));
        else
            dr_sub = 0;
        end
        
        % Sub-pixel offset in columns
        if x(1) ~= x(3)
            dc_sub = (x(3) - x(1)) / (2 * (2 * x(2) - x(1) - x(3)));
        else
            dc_sub = 0;
        end

        % Convert to signed displacement including sub-pixel component
        dR = (rp - 1 + dr_sub); if dR > rows/2, dR = dR - rows; end
        dC = (cp - 1 + dc_sub); if dC > cols/2, dC = dC - cols; end
        
        shifts(i, :) = [dC, dR];
    end
    
    % 3. Apply Correction & Auto-Crop
    % Calculate max jitter to hide black borders
    maxExcursion = ceil(max(abs(shifts(:))));
    cropMargin = maxExcursion + 5; 
    rect = [cropMargin, cropMargin, cols - 2*cropMargin, rows - 2*cropMargin];
    
    % v.CurrentTime = 0;
    % tempFrame = readFrame(v);
    tempFrame = refFrame;
    tempCorrected = imtranslate(tempFrame, shifts(1, :), 'linear', 'OutputView', 'same');
    tempFinal = imcrop(tempCorrected, rect);
    [finalH, finalW, channels] = size(tempFinal);
    correctedFrames = zeros(finalH, finalW, channels, totalFrames, 'uint8');

    fprintf('Applying stabilization...\n');
    
    for i = 1:totalFrames
        % if ~hasFrame(v), break; end
        frame = refFrame;
        
        % Translate and Crop
        % Sub-pixel precision is handled by imtranslate
        corrected = imtranslate(frame, shifts(i, :), 'linear', 'OutputView', 'same');
        finalFrame = imcrop(corrected, rect);
        correctedFrames(:,:,:,i) = finalFrame;
        % writeVideo(vOut, finalFrame);
    end
    
    % close(vOut);
    disp('Done! Video is now rock-solid.');
end
