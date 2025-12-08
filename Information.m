%% transducer parameters; 
% For more info: 
% https://cw.fel.cvut.cz/b232/_detail/courses/zsl/schematization-of-a-linear-array-with-the-width-height-pitch-and-kerf-parameters-of.png?id=courses%3Abam33zsl%3Alabs2023_08_ussim

data = load("RcvData.mat");
rf = data.RcvData;

probe_central_frequency = 2.841*1e6;% unit Hz.
number_of_elements = 80;            % number of elements on the probe.
element_width = 2.3*1e-4;           % unit m; width of each individual element.
element_pitch = 2.7*1e-4;           % unit m; pitch of each individual element.
element_kerf = 4.0*1e-5;            % unit m; kerf of each individual element.
probe_band_width =[1500000,4200000];% unit Hz; the bandwidth of the probe.

%% wave transmit parameters
transmit_frequency = 2.0*1e6;       % unit Hz; this is the frequency we used for imaging, it can be different from the probe's central frequency.
imaging_depth = 80*1e-3;            % unit m; the deepest position of the recorded signal.
sampling_rate = 14.205*1e6;         % unit Hz; the sampling rate of the recorded signal.
steering_angle = 0;                 % we used only one plane wave for imaging. No steered plane wave was used.
speed_of_sound = 1540;              % unit m; the speed of sound.

%% image reconstruction parameters
pixelmap_x = -0.06:(2.7*1e-4):0.06; % unit m, the 'x' coordinates of of reconstruction area;
pixelmap_z = 0:(3.4*1e-4):0.08;     % unit m, the 'z' coordinates of of reconstruction area;
% we are doing 2D imaging, thus there is no 'y' coordinates, or y is equals to 0.


%% dataset infomation
% the size of the RcvData.mat file RcvData is 1540 * 80.
% 1540 is total number of samples in each channel (or element).
% 80 is the total number of elements.
% Each column of the array indicates a *RF signal* recorded from the elements. (not IQ, you will read it in the DAS paper).
% You may find out the the total number of samples matches the following
% calculation: 2*hypot(transducer_total_length, imaging_depth)/speed_of_sound*sampling_rate.

% the above info should be enough for you to reconstruct one image frame.
% if you need more information, please let us know.
