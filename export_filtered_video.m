function export_filtered_video(db_stack, outputFile, displayFloor, frameRate)
% db_stack: HWT dB stack
% outputFile: output mp4 filename
% displayFloor:-35 or -30
% frameRate: 15

    if nargin < 4
        frameRate = 15;
    end

    v = VideoWriter(outputFile, 'MPEG-4');
    v.FrameRate = frameRate;
    open(v);

    fig = figure('Position',[100 100 450 850], 'Color','black');

    for t = 1:size(db_stack,3)
        clf
        img = db_stack(:,:,t);
        img(img < displayFloor) = displayFloor;

        imagesc(img, [displayFloor 0]);
        axis image off
        colormap gray
        set(gca, 'Color', 'black')
        drawnow

        fr = getframe(fig);
        writeVideo(v, fr);
    end

    close(v);
end