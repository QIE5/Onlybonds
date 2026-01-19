%% Discretation calculations
% ≥0.5λ in the lateral direction , ≥ 0.25λ in the axial direction.
fc = 2.0*1e6;
c = 1540;
lambda = c /fc; % λ = c / f
pixelSize = 3.3879*1e-5;
lambda/4 > pixelSize; % This returns true, so discretation is enough



%% Code for producing video with localied MBs
simvid = VideoReader('simulation.mp4');
% localisedvid = VideoWriter("localised");
% localisedvid.FrameRate = simvid.FrameRate;
% open(localisedvid);
% figure('Visible','off');

numFrames = simvid.NumFrames;
centroidCell = cell(numFrames, 1);
frame1 = readFrame(simvid);
figure;
imshow(frame1);
hold on
plot(417.4997 ,1.2695e+03,'b*');
% The lines below are drawn through the centre of the selected MB.
% They represent the dimensions of the MB
plot ([417.4997 - 90 ,417.4997 + 90], [ 1.2695e+03,1.2695e+03], 'r-'); % PSF width = 180
plot ([417.4997,417.4997], [ 1.2695e+03 + 15,1.2695e+03 - 15], 'g-'); % heaight = 30
hold off
patchSize = 180 + 3;
%%
numFrames = 2;
for n = 1: numFrames;
    frame = read(simvid, n);
    frame = im2gray(frame);
    threshold = prctile(frame(:), 99.5);  % top 0.5% intensities, arbitrarily selected to include starburst
    regionalMax = imregionalmax(frame);
    % bw = frame > threshold;
    bw = regionalMax & (frame > threshold);
    % [xPeaks, yPeaks] = find(bw)
    
    cc = bwconncomp(bw, 8);
    stats = regionprops(cc, frame, 'WeightedCentroid');

    centroidCell{n} = cat(1,stats.WeightedCentroid);

    imshow(frame); hold on
    plot(centroidCell{n}(:,1), centroidCell{n}(:,2), 'r+', 'MarkerSize',6)
    hold off

    % F = getframe(gca);
    % writeVideo(localisedvid, F);
end
%%
centroids = centroidCell{1};
figure;
imshow(frame1);
hold on
plot(centroids(:,1),centroids(:,2),'b*');
hold off
title('Frame 1');

% close(localisedvid);
% close(gcf);



%lsqcurvefit
