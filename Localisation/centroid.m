simvid = VideoReader('simulation.mp4');


frame = readFrame(simvid);
frame = im2gray(frame);
threshold = prctile(frame(:), 95)  % top 0.5% intensities, arbitrarily selected to include starburst
bw = frame > threshold;
cc = bwconncomp(bw, 8);
stats = regionprops(cc, frame, 'WeightedCentroid');
centroids = cat(1,stats.WeightedCentroid);


figure;
imshow(frame);
hold on
plot(centroids(:,1),centroids(:,2),'b*')
hold off
title('Frame 1');
