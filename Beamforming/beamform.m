function [bfImageDB] = beamform(rawSIG, beamformParam, varargin)
% beamform - Simple DAS wrapper that handles RF or IQ input
%
% Usage:
%   bfImageDB = beamform(rawSIG, beamformParam)
%   bfImageDB = beamform(rawSIG, beamformParam, 'SignalType', 'RF');

    %% 1. Input parser
    defaultSignalType = 'IQ'; 
    validSignalTypes = {'RF','IQ'};
    
    p = inputParser;
    addParameter(p, 'SignalType', defaultSignalType, @(x) any(validatestring(x, validSignalTypes)));
    parse(p, varargin{:});                   % <-- must parse before accessing results
    signalType = p.Results.SignalType;

    %% 2. Convert RF → IQ if needed
    if strcmp(signalType, 'RF')
        iq = rf2iq(double(rawSIG), beamformParam);
    else
        iq = rawSIG;
    end


    % Using IQ signal (MUST toolbox):
    pixelmapX = beamformParam.pixelmapX;
    pixelmapZ = beamformParam.pixelmapZ;
    [X, Z] = meshgrid(pixelmapX, pixelmapZ);
    vsource = [0, -10000]; % default, can add functionality later
    % x0=steering_angle=0, z0 is large and negative 
    % if there is no steering angle, then the virtual source is directly above
    % the array. This makes theta perpendicular to the z axis
    % beta can also be calculated based on the virtual source
    [bfSIG, M] = ezdas(iq, X(:), Z(:), vsource, beamformParam);
    
    % Using Rf signal:
    % [bfSIG, M] = ezdas(double(rf), X(:), Z(:), vsource, param);
    
    % Reshape to image dimensions
    bfImage = reshape(bfSIG, size(X));
    
    % Find the envelope
    
    % Using IQ signal:
    envelope = abs((bfImage));
    
    % Using Rf signal:
    % envelope = abs(hilbert((bfImage)));

    bfImageDB = 20*log10(envelope + eps);
    bfImageDB = bfImageDB - max(bfImageDB(:));
    bfInt = mat2gray(bfImageDB, [-40 0]) * 255;
end

% This uses ezdas.m and the MUST toolbox



%% Display the result
% figure('Position', [100, 100, 800, 600]);
% 
% % Display in dB scale
% 
% imagesc(pixelmapX, pixelmapZ, bfImageDB, [-40, 0]); 
% colormap('gray');
% colorbar;
% xlabel('Lateral Distance (m)');
% ylabel('Axial Distance (m)');
% title('DAS Beamformed Image');
% axis image;

