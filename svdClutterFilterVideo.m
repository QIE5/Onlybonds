function filteredStack = svdClutterFilterVideo(videoInput, nCut, nNoiseCut)
%Supported inputs:VideoReader object+ HWT RGB frame stack+ HWT grayscale frame stack
%ideal filteredStack is HWT filtered grayscale stack

    if isa(videoInput, 'VideoReader')
        frames = read(videoInput);

        if ndims(frames) ~= 4 || size(frames, 3) ~= 3
            error('VideoReader input must produce RGB frames of size H x W x 3 x T');
        end

        T = size(frames, 4);
        grayStack = zeros(size(frames,1), size(frames,2), T, 'double');

        for t = 1:T
            grayFrame = rgb2gray(frames(:,:,:,t));
            grayStack(:,:,t) = im2double(grayFrame);
        end

    else
        if ~isnumeric(videoInput)
            error('Input must be either a VideoReader object or a numeric frame stack');
        end

        if ndims(videoInput) == 4
            % situation2
            if size(videoInput,3) ~= 3
                error('4D numeric input must be HxWx3xT');
            end

            T = size(videoInput, 4);
            grayStack = zeros(size(videoInput,1), size(videoInput,2), T, 'double');

            for t = 1:T
                grayFrame = rgb2gray(videoInput(:,:,:,t));
                grayStack(:,:,t) = im2double(grayFrame);
            end

        elseif ndims(videoInput) == 3
            % situation3
            grayStack = im2double(videoInput);

        else
            error('Numeric input must be HxWxT or HxWx3xT');
        end
    end

    filteredStack = svdClutterFilter(grayStack, nCut, nNoiseCut);
end