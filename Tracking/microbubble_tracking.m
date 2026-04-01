function tracks = microbubble_tracking(detections, dt, params)
%   detections: 1xT cell, detections{t} is Nt x 2 matrix [x z].
%   dt: time between frames (seconds).
%   params: struct with tracking parameters（as containers)
%   tracks: struct array of finished tracks.

     %0) Basic input validation

    if nargin < 2 % check function have enough parameters
        error('not enough inpt parameters');
    end
    if nargin < 3 && isempty(params) %parameters empty then consturct one array
        params = struct();
    end
    if ~iscell(detections) %check detections is a cell array or error
        error('detections must be a cell array like detections{t} = [x z].');
    end
    T = numel(detections); %count how many frames in total
    if T == 0
        tracks = struct([]);
        return;
    end


    % 1) Set defaults 
%set params as a box to store the data
    if ~isfield(params,'gate')
        params.gate = 12;  % max allowed displacement per frame (same units as x/z)
    end
    if ~isfield(params,'max_missed')
        params.max_missed = 2; % allow track to be missing this many frames
    end
    if ~isfield(params,'min_length')
        params.min_length = 3; % filter out very short tracks at the end
    end
    if ~isfield(params,'use_prediction')
        params.use_prediction = true; %constant-velocity prediction
    end
    if ~isfield(params,'vel_smooth')
        params.vel_smooth = 0.8; %larger = smoother velocity
    end
    if ~isfield(params,'costOfNonAssignment')
        params.costOfNonAssignment = (params.gate^2); % penalty for leaving a track/detection unmatched
    end
    if ~isfield(params,'bigM')
    params.bigM = 1e9;
    end
   
    % 2) Track state containers
 

    % 2) Track state containers (safe empty struct with fields)
    activeTracks   = struct('id',{}, 'frames',{}, 'x',{}, 'z',{},'vx',{}, 'vz',{}, 'hasVelocity',{}, 'missed',{});
    finishedTracks = struct('id',{}, 'frames',{}, 'x',{}, 'z',{},'vx',{}, 'vz',{}, 'hasVelocity',{}, 'missed',{});
    nextId = 1;

    %activeTracks = struct([]);% create activetracks to store exist microbbules/initialise
    %finishedTracks = struct([]);%create finishedtracks to store microbubbles that leave the detecting region
    %nextId = 1; %create nextid to give each tracks unique id

    % 3) Main loop over frames

    for t = 1:T
        Z = detections{t}; %extract the cell on number t from T
        if isempty(Z) %frames that don't have microbubbles
            Z = zeros(0,2);%to create empty matrix even this frame don't have microbubble to prevent system breaking
        end
        if size(Z,2) ~= 2 % to make sure Z have 2 columns
            error('detections{%d} must be Nt x 2: [x z].', t);
        end
        nDet = size(Z,1);%how many microbubbles detect in this frame
        nTrk = numel(activeTracks);%to count how many activetracks

       
        % 3.1) If no active tracks, start new tracks from detections
        
        if nTrk == 0
            for j = 1:nDet %go through all the microbubbles
                activeTracks(end+1) = newTrack(nextId, t, Z(j,1), Z(j,2)); %Z(j,1)for j microbubbles xcoordinate, Z(j,2) for j microbubbles ycoordinate, t for frames， initialize
                nextId = nextId + 1;
            end
            continue;%to execute for t=1:T
        end

        % 3.2) Predict track positions at frame t

        pred = zeros(nTrk,2);%create a new matrix with dimension nTrk*2
        for i = 1:nTrk % to make prediction for every activetracks
            lastX = activeTracks(i).x(end);% The x-coordinate of the i-th trajectory last observed
            lastZ = activeTracks(i).z(end);%same for z-coordinate
            if params.use_prediction && activeTracks(i).hasVelocity %to ensure that not the first microbubble case,appear twice can have velocity
                pred(i,1) = lastX + activeTracks(i).vx * dt; %predict xcoordinate
                pred(i,2) = lastZ + activeTracks(i).vz * dt; %predict zcoordinate
            else
                pred(i,1) = lastX;
                pred(i,2) = lastZ;
            end
        end

        % 3.3) Build cost matrix (tracks x detections), apply gating
        gate2 = params.gate^2;

        C = zeros(nTrk, nDet); %to create matrix where tracks*points and (ith,jth)means ith tracks jth microbubbles
        for i = 1:nTrk %go through all the tracks
            dx = Z(:,1) - pred(i,1); %difference between real detection xcoordinate and prediction one
            dz = Z(:,2) - pred(i,2);
            C(i,:) = (dx.^2 + dz.^2).';
        end

        C(C > gate2) = params.bigM; %has two way, first is the bigM, second is costOfNonAssignment=gate^2 to solve match or unmatch problems 

   
        
        
     % 3.4) Assignment (Hungarian)
        %     using assignDetectionsToTracks(costMatrix, costOfNonAssignment)
        if nDet == 0  %no microbubbles detected
            assignments = zeros(0,2); %assignment is [trackIndex,detectionIndex]
            unassignedTracks = (1:nTrk).'; % to list all the tracks that is not matched with microbubbles
            unassignedDetections = zeros(0,1); % to list all the microbubbles that is not matched with tracks
        else
            [assignments, unassignedTracks, unassignedDetections] = assignDetectionsToTracks(C, params.costOfNonAssignment);
        end

        % 3.5) Update matched tracks
        for k = 1:size(assignments,1) %how many rows assignments have, that is successful matched pairs (detection and track)
            trkIdx = assignments(k,1); %k-th row 1st column
            detIdx = assignments(k,2); %k-th row 2nd column

            xNew = Z(detIdx,1); % x-coordinate of microbubbles
            zNew = Z(detIdx,2); % z-coordinate of microbubbles

            xPrev = activeTracks(trkIdx).x(end); %tracks last observed location to calculate velocity
            zPrev = activeTracks(trkIdx).z(end);

            activeTracks(trkIdx).frames(end+1) = t; %append new frame to activetracks
            activeTracks(trkIdx).x(end+1) = xNew; %append xnew coordinate to the end of tracks x 
            activeTracks(trkIdx).z(end+1) = zNew; %same for z

            vxNew = (xNew - xPrev) / dt;
            vzNew = (zNew - zPrev) / dt;
            if activeTracks(trkIdx).hasVelocity
                a = params.vel_smooth; %give value
                activeTracks(trkIdx).vx = a*activeTracks(trkIdx).vx + (1-a)*vxNew;
                activeTracks(trkIdx).vz = a*activeTracks(trkIdx).vz + (1-a)*vzNew;
            else
                activeTracks(trkIdx).vx = vxNew;
                activeTracks(trkIdx).vz = vzNew;
                activeTracks(trkIdx).hasVelocity = true;
            end

            activeTracks(trkIdx).missed = 0; %successfully match microbubbles at this frame so missed START from 0
        end


        % 3.6) Handle unpaired tracks (missed detections)

        toFinish = false(nTrk,1); %create all false list/vector
        for ii = 1:numel(unassignedTracks)
            trkIdx = unassignedTracks(ii);
            activeTracks(trkIdx).missed = activeTracks(trkIdx).missed + 1;
            if activeTracks(trkIdx).missed > params.max_missed
                toFinish(trkIdx) = true;
            end
        end

        if any(toFinish)
            finishedTracks = [finishedTracks, activeTracks(toFinish)]; %all missed tracks come together 
            activeTracks = activeTracks(~toFinish); %extract tracks that not missed, (still active)
        end


        % 3.7) Create new tracks for unassigned detections (new microbubbles)

        for ii = 1:numel(unassignedDetections)
            detIdx = unassignedDetections(ii);
            activeTracks(end+1) = newTrack(nextId, t, Z(detIdx,1), Z(detIdx,2)); %same format as newTrack above, to handle new birth microbubbles
            nextId = nextId + 1;
        end
    end

    % 4) Finalize: move all remaining active to finished
    finishedTracks = [finishedTracks, activeTracks];

    % 5) Filter tracks by minimum length
    keep = false(1, numel(finishedTracks));
    for i = 1:numel(finishedTracks)
        keep(i) = numel(finishedTracks(i).frames) >= params.min_length; % to remove the noise that only exist min_length frames
    end
    tracks = finishedTracks(keep);
end


% function that create new track

function trk = newTrack(id, t, x, z)
    trk.id = id;
    trk.frames = t;
    trk.x = x;
    trk.z = z;
    trk.vx = 0;
    trk.vz = 0;
    trk.hasVelocity = false;
    trk.missed = 0;
end