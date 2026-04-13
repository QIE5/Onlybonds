% Each extracted PSF is stored here. To use, copy the respective block into
% psfFinder.m


%% Simulation.mp4 microbubbles
simvid = VideoReader('simulation.mp4');
refFrame1 = read(simvid, 1);
refFrame2 = read(simvid, 30);
figure;
imshow(refFrame1);
impixelinfo;
hold on
bubble1 = [508,1068]; % Frame 1
plot(bubble1(1), bubble1(2), 'b*');
bubble2 = [435,1489]; % Frame 30
% plot(bubble2(1), bubble2(2), 'b*');
bubble3 = [534, 2232]; % Frame 1
plot(bubble3(1), bubble3(2), 'b*');

% The lines below are drawn through the centre of the selected MB.
% They represent the dimensions of the MB
% Bubble 1
plot ([bubble1(1) - 100 ,bubble1(1) + 100], [bubble1(2), bubble1(2)], 'r-'); % PSF width = 200
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + 15,bubble1(2) - 15], 'g-'); % height = 30
psfWidth1 = 200;
psfHeight1 = 30;
% Bubble 2, frame 30
plot ([bubble2(1) - 85 ,bubble2(1) + 85], [bubble2(2), bubble2(2)], 'r-'); % PSF width = 170
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 14,bubble2(2) - 14], 'g-'); % height = 28
psfWidth2 = 170;
psfHeight2 = 28;
% Bubble 3
plot ([bubble3(1) - 70 ,bubble3(1) + 40], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 110
% plot ([bubble3(1), bubble3(1)], [ bubble3(2) + 12,bubble3(2) - 12], 'g-'); % height = 24
psfWidth3 = 110;
psfHeight3 = 24;


% Clutter filter video
simvid = VideoReader('static_background_clutter_filterd.mp4');
refFrame1 = read(simvid, 70);
refFrame2 = read(simvid, 140);
figure;
imshow(refFrame1);
title ('Display PSF regions')
hold on

% Bubble 1
bubble1 = [113,44]; % Frame 70
plot(bubble1(1), bubble1(2), 'b*');
plot ([bubble1(1) - 35 ,bubble1(1) + 35], [bubble1(2), bubble1(2)], 'r-'); % PSF width = 70
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + 6,bubble1(2) - 6], 'g-'); % height = 12
psfWidth1 = 70;
psfHeight1 = 12;

% Bubble 2
bubble2 = [94,151]; % Frame 70
plot(bubble2(1), bubble2(2), 'b*');
% Bubble 2, frame 140
plot ([bubble2(1) - 34 ,bubble2(1) + 34], [bubble2(2), bubble2(2)], 'r-'); % PSF width = 64
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 6,bubble2(2) - 6], 'g-'); % height = 12
psfWidth2 = 64;
psfHeight2 = 12;

% Bubble 3
bubble3 = [82,244]; % Frame 70
plot(bubble3(1), bubble3(2), 'b*');
plot ([bubble3(1) - 34 ,bubble3(1) + 34], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 68
plot ([bubble3(1), bubble3(1)], [ bubble3(2) + 6,bubble3(2) - 5], 'g-'); % height = 11
psfWidth3 = 68;
psfHeight3 = 11;

%% Clutter filter video
simvid = VideoReader('static_background_clutter_filterd.mp4');
refFrame1 = read(simvid, 70);
refFrame2 = read(simvid, 140);
figure;
imshow(refFrame1);
impixelinfo;
hold on

% Bubble 1
bubble1 = [113,44]; % Frame 70
plot(bubble1(1), bubble1(2), 'b*');
plot ([bubble1(1) - 35 ,bubble1(1) + 35], [bubble1(2), bubble1(2)], 'r-'); % PSF width = 70
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + 6,bubble1(2) - 6], 'g-'); % height = 12
psfWidth1 = 70;
psfHeight1 = 12;

% Bubble 2
bubble2 = [94,151]; % Frame 70
plot(bubble2(1), bubble2(2), 'b*');
% Bubble 2, frame 140
plot ([bubble2(1) - 34 ,bubble2(1) + 34], [bubble2(2), bubble2(2)], 'r-'); % PSF width = 64
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + 6,bubble2(2) - 6], 'g-'); % height = 12
psfWidth2 = 64;
psfHeight2 = 12;

% Bubble 3
bubble3 = [82,244]; % Frame 70
plot(bubble3(1), bubble3(2), 'b*');
plot ([bubble3(1) - 34 ,bubble3(1) + 34], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 68
plot ([bubble3(1), bubble3(1)], [ bubble3(2) + 6,bubble3(2) - 5], 'g-'); % height = 11
psfWidth3 = 68;
psfHeight3 = 11;

%% CEUS Stable2
simvid = VideoReader('Phantom Videos/CEUS_Stable2.mp4');
numframes = simvid.NumFrames;
refFrame1 = read(simvid, 810);
refFrame1 = thresholding(refFrame1);
refFrame1 = im2gray(refFrame1);
refFrame2 = read(simvid, 270);
refFrame2 = thresholding(refFrame2);
refFrame3 = read(simvid, 390);
refFrame3 = thresholding(refFrame3);

displayFrame = refFrame3; % Pick the display frame
figure;
imshow(displayFrame);
impixelinfo;
zoom on;
title ('Display PSF regions')
hold on

% Bubble 3
bubble3 = [1448,589]; % Frame 810, refFrame1
plot(bubble3(1), bubble3(2), 'b*');
plot ([bubble3(1) - 32 ,bubble3(1) + 38], [bubble3(2), bubble3(2)], 'r-'); % PSF width = 70
plot ([bubble3(1), bubble3(1)], [ bubble3(2) + 7,bubble3(2) - 7], 'g-'); % height = 14
psfWidth1 = 70;
psfHeight1 = 14;

% Bubble 2
bubble2 = [662,604]; % refFrame2
plot(bubble2(1), bubble2(2), 'b*');
psfWidth2 = 70;
psfHeight2 = 24;
plot ([bubble2(1) - psfWidth2/2 ,bubble2(1) + psfWidth2/2], [bubble2(2), bubble2(2)], 'r-');%
plot ([bubble2(1), bubble2(1)], [ bubble2(2) + psfHeight2/2,bubble2(2) - psfHeight2/2], 'g-'); 


% Bubble 1
bubble1 = [362,605]; % refFrame3
plot(bubble1(1), bubble1(2), 'b*');
psfWidth3 = 60;
psfHeight3 = 18;
plot ([bubble1(1) - psfWidth3/2 ,bubble1(1) + psfWidth3/2], [bubble1(2), bubble1(2)], 'r-'); 
plot ([bubble1(1), bubble1(1)], [ bubble1(2) + psfHeight3/2,bubble1(2) - psfHeight3/2], 'g-'); 