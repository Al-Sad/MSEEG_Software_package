%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Copyright 2019 Mohammad Al-Sa'd
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Mohammad F. Al-Sa'd     (mohammad.al-sad@tuni.fi)
%          Prof. Boualem Boashash  (b.boashash@uq.edu.au)
%
% The following reference should be cited whenever this script is used:
%   M. Al-Sa'd, B. Boashash, "Design and implementation of a multi-sensor
%   newborn EEG seizure and background model with inter-channel field
%   characterization", Digital Signal Processing, 2019.
%
% Last Modification: 25-Jan-2019
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%            10-20 International standard electrode placement
%
%  Syntax  : [p,tag] = electrode_10_20(R)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R        : Head radius in cm.
%
% <OUTPUTs>
% p        : Cartesian 3D position of each EEG electrode.
% tag      : Label of each EEG electrode.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% [p,tag] = electrode_10_20(5);
% disp({'Label' 'x' 'y' 'z'})
% disp([char(tag) num2str(p)]);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [p,tag] = electrode_10_20(R)
%% NOTE
% Each electrode has its position initially in spherical coordinates
% as follows: [Azimuth, Elevation] 
% Azimuth angle is measured from the +X
% Elevation angle is measured from the X-Y Projection to the Vector

%% Electrode Labelling
tag = {'F_z';'C_z';'P_z';'F_p_1';'F_3';'F_7';'C_3';'T_3';'P_3';'T_5';'O_1';...
       'F_p_2';'F_4';'F_8';'C_4';'T_4';'P_4';'T_6';'O_2';'F_p_z';'O_z'};
%% Center Electrodes
Fz = [180 ; 54]; Cz = [180 ; 90]; Pz = [360 ; 54];
center_elec = [Fz Cz Pz]*pi/180;
[p(1:3,1), p(1:3,2), p(1:3,3)] = sph2cart(center_elec(1,:), center_elec(2,:),R);

%% Left-Side Electrodes
Fp1 = [198 ; 18];      F3 = [218.36 ; 42.3]; F7 = [234 ; 18];
C3  = [270 ; 54];      T3 = [270 ; 18];
P3  = [321.64 ; 42.3]; T5 = [306 ; 18];      O1 = [342 ; 18];
left_elec = [Fp1 F3 F7 C3 T3 P3 T5 O1]*pi/180;
[p(4:11,1), p(4:11,2), p(4:11,3)] = sph2cart(left_elec(1,:), left_elec(2,:),R);

%% Right-Side Electrodes
Fp2 = [162 ; 18];     F4 = [141.64 ; 42.3]; F8 = [126 ; 18];
C4  = [90 ; 54];      T4 = [90 ; 18];
P4  = [38.36 ; 42.3]; T6 = [54 ; 18];       O2 = [18 ; 18];
right_elec = [Fp2 F4 F8 C4 T4 P4 T6 O2]*pi/180;
[p(12:19,1), p(12:19,2), p(12:19,3)] = sph2cart(right_elec(1,:), right_elec(2,:),R);

%% Front and Back Electrodes
Fpz = [180 ; 18]; Oz = [360 ; 18];
front_back_elec = ([Fpz Oz]*pi/180);
[p(20:21,1), p(20:21,2), p(20:21,3)] = sph2cart(front_back_elec(1,:), front_back_elec(2,:),R);

end
