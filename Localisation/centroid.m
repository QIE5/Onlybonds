simvid = VideoReader('simulation.mp4');
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
end
centroids = centroidCell{1};
figure;
imshow(frame1);
hold on
plot(centroids(:,1),centroids(:,2),'b*');
hold off
title('Frame 1');
