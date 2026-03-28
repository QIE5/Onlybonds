function [stack_static, db_static] = run_static_filter(staticframe)

%staticframe: HWT beamformed data
%stack_static: filtered complex stack
%db_static: normalized dB display stack

    nCut = 4;
    nNoiseCut = 0;

    stack_static = svdClutterFilter(staticframe, nCut, nNoiseCut);

    env_static = abs(stack_static);
    env_static = env_static ./ max(env_static(:));
    db_static = 20*log10(env_static + eps);

end