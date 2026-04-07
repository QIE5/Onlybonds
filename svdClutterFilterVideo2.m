function filteredFrames = svdClutterFilterVideo2(frames, nCut, nNoiseCut)
% Input:
%   frames: H x W x T grayscale video stack
% Output:
%   filteredFrames: H x W x T filtered video stack

    frames = im2double(frames);
    stack_filt = svdClutterFilter(frames, nCut, nNoiseCut);
    filteredFrames = abs(stack_filt);
end