function [localisedBubbleCoords] = crossCorrelation(frame, localisationParam)
% function [localisedBubbleCoords] = crossCorrelationLocal(frame, localisationParam)
    [height, width, ~] = size(frame);
    viableYTop = 570; 
    viableYBottom = 850;
    % viableYTop = 1; Simulation
    % viableYBottom = height;
    % Hard-coded, where the bubbles start on the 
    % frame. Ideally, we should only be passing the relevant part of the 
    % frame into the function. Minimum value is 1.
    
    %% Thresholding
    threshold = prctile(frame(:), 95); % Adding thresholding to function
    frame(frame <= threshold) = 0;
    %% Detection
    
    psfTemplates = localisationParam.psfTemplates;
    
    % edges = round(linspace(viableYTop, height, 4)); % simulation
    edges = round(linspace(1, width, 4));
    allCoords = cell(numel(psfTemplates),1);
    roiYMin = viableYTop;
    roiYMax = viableYBottom;
    for n = 1:numel(psfTemplates)
        % roiYMin = edges(n);
        % roiYMax = edges(n+1); % Uncomment for simulation
        % roi = frame(roiYMin:roiYMax, :);

        roiXMin = edges(n);
        roiXMax = edges(n+1);
        roi = frame(roiYMin: roiYMax, roiXMin:roiXMax);
   

        c = normxcorr2(psfTemplates{n}, roi);
    
        th = round(size(psfTemplates{n},1)/2);
        tw = round(size(psfTemplates{n},2)/2);
        
        c_valid = c(th:end-th+1, tw:end-tw);
        % corc = c_valid > 0.7; % Simulation
        corc = c_valid > 0.7;
        shading flat
        stats = regionprops(corc, c_valid, 'Area', 'WeightedCentroid', 'BoundingBox');
        if ~isempty(stats)
            coords = cat(1, stats.WeightedCentroid);
            coords(:,2) = coords(:,2) + roiYMin - 1;
            coords(:,1) = coords(:,1) + roiXMin - 1;
            % coords(:,1) = coords(:,1);
            allCoords{n} = coords;
        end
    end
    localisedBubbleCoords = vertcat(allCoords{:});
    
end

% Uncomment to run code locally
% bubbleVid = VideoReader('simulation.mp4');
% bubbleVid = VideoReader('static_background_clutter_filterd.mp4');
% bubbleVid = VideoReader('Phantom Videos/CEUS_Stable1.mp4');
% 
% numFrames = bubbleVid.NumFrames;
% frame = read(bubbleVid, 390);
% frame = im2gray(frame);
% 
% localisationParam.psfTemplates = {psfTemplate1, psfTemplate2, psfTemplate3};
% [localisedBubbleCoords] = crossCorrelationLocal(frame, localisationParam);
% 
% % View singular frame
% threshold = prctile(frame(:), 90); % Adding thresholding to function
% frame(frame <= threshold) = 0;
% figure;
% imshow(frame);
% impixelinfo;
% hold on
% plot(localisedBubbleCoords(:,1),localisedBubbleCoords(:,2),'b*');
% hold off
% title('Frame 1');

% %% Reconstruct video
% localisedvid = VideoWriter("stable2");
% localisedvid.FrameRate = bubbleVid.FrameRate;
% open(localisedvid);
% fig = figure('Color','k');
% ax  = axes(fig);
% ax.Position = [0 0 1 1]; 
% localisedBubbleCoords = cell(numFrames, 1);
% 
% for n = 1:numFrames
%     frame = read(bubbleVid, n);
%     frame = im2gray(frame);
%     [localisedBubbleCoords{n}] = crossCorrelationLocal(frame, localisationParam);
% 
%     imshow(frame, 'Parent', ax); hold(ax, 'on');
%     plot(ax, localisedBubbleCoords{n}(:,1), localisedBubbleCoords{n}(:,2), ...
%         'b*');
%     hold(ax, 'off');
% 
%     F = getframe(ax);          % capture AXES only, not figure
%     writeVideo(localisedvid, F);
% end
% 
% 
% close(localisedvid);