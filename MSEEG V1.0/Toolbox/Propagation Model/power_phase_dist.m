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
%                       Power-Phase distributions
%
% Syntax  : [P, PHI, p] = power_phase_dist(event_p, N, fs, mP, v)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% event_p : EEG source signals positions in 3D cartesian coordinates within
%           the brain sphere. It is a Px3 matrix, where P is the number of
%           EEG sources.
% N       : Total number of samples for the distribution.
% fs      : Sampling frequency in Hz.
% mP      : Mean multi-sensor power of EEG.
% v       : Propagation speed of source signals in cm/s.
%
% <OUTPUTs>
% P       : Power distribution on the scalp.
% PHI     : Phase distribution on the scalp.
% p       : Scalp manifold.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% event_p = [-2.7 2.7 1.2 ; 0 -3 2 ; 4 0 0];
% N = 256; fs = 32; mP = 100; v = 5;
% [P, A, p] = power_phase_dist(event_p, N, fs, mP, v);
% figure; colormap jet; title('Power distribution');
% surface(p{1}, p{2}, p{3}, P,'edgealpha',0);
% figure; colormap jet; title('Phase distribution');
% surface(p{1}, p{2}, p{3}, A,'edgealpha',0);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [P, PHI, p] = power_phase_dist(event_p, N, fs, mP, v)
%% Parameters
Scalp_R  = 5.95;        % Head Radius in cm

%% Computing Scalp locations
azi = linspace(0, 2*pi, N);
ele = linspace(0, pi/2, N);
[Azi, Ele] = meshgrid(azi, ele);
[p{1},p{2},p{3}] = sph2cart(Azi, Ele, Scalp_R);

%% Event Generation and amplitude Dispersion
V = zeros(N,N);
PHI = zeros(N,N);
for i = 1:N
    for j = 1:N
        [V(i,j), PHI(i,j),] = multi_sensor_propagation(event_p, fs, mP, v, [p{1}(i,j), p{2}(i,j), p{3}(i,j)]);
    end
end

%% Power Dispersion
P = V.^2;
PHI = PHI./fs;
end