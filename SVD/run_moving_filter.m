function [stack_moving, db_moving, int_moving] = run_moving_filter(movingframe)

%movingframe: HWT beamformed data
%stack_moving: filtered complex stack
%db_moving: normalized dB display stack

    nCut = 5;
    nNoiseCut = 0;

    stack_moving = svdClutterFilter(movingframe, nCut, nNoiseCut);

    env_moving = abs(stack_moving);
    env_moving = env_moving ./ max(env_moving(:));
    db_moving = 20*log10(env_moving + eps);
    int_moving = mat2gray(db_moving, [-40 0]) * 255;
end