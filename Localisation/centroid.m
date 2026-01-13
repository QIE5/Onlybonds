simvid = VideoReader('simulation.mp4');

%% Code for producing video with localied MBs
% localisedvid = VideoWriter("localised");
% localisedvid.FrameRate = simvid.FrameRate;
% open(localisedvid);
% figure('Visible','off');

numFrames = simvid.NumFrames;
centroidCell = cell(numFrames, 1);
frame1 = readFrame(simvid);


for n = 1: numFrames;
    frame = read(simvid, n);
    frame = im2gray(frame);
    threshold = prctile(frame(:), 99.5);  % top 0.5% intensities, arbitrarily selected to include starburst
    bw = frame > threshold;
    cc = bwconncomp(bw, 8);
    stats = regionprops(cc, frame, 'WeightedCentroid');
    centroidCell{n} = cat(1,stats.WeightedCentroid);
    %%
    % imshow(frame); hold on
    % plot(centroidCell{n}(:,1), centroidCell{n}(:,2), 'r+', 'MarkerSize',6)
    % hold off
    % 
    % F = getframe(gca);
    % writeVideo(localisedvid, F);
end
centroids = centroidCell{1};
figure;
imshow(frame1);
hold on
plot(centroids(:,1),centroids(:,2),'b*');
hold off
title('Frame 1');

% close(localisedvid);
% close(gcf);