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
%                           Ears polygon model
%
%  Syntax  : [p_r,p_l] = ellipsoid_model(R, azi_s, azi_f, N, M)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R        : Head radius in cm to center the two ellipsoids.
% azi_s    : Start of azimuth spanning angle.
% azi_f    : End of azimuth spanning angle.
% N        : Meshing resolution of the elevation angle.
% M        : Meshing resolution of the azimuth angle.
%
% <OUTPUTs>
% p_r      : Right ear polygon in cartesian coordinates.
% p_l      : Left ear polygon in cartesian coordinates.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% [p_r,p_l] = ellipsoid_model(4, 0, 2*pi, 32, 32);
% surf(p_r{1},p_r{2},p_r{3}); hold on;
% surf(p_l{1},p_l{2},p_l{3});
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [p_r, p_l] = ellipsoid_model(R, azi_s, azi_f, N, M)
azi = linspace(azi_s, azi_f, M);
ele = linspace(-pi/2, pi/2, N);
[Azi, Ele] = meshgrid(azi, ele);

%% Right Ear
R_x = 1; R_y = 0.5; R_z = 0.5;
[ear_r_x, ear_r_y, ear_r_z] = sph2cart(Azi, Ele, 1);
p_r{1} = R_x.*ear_r_x;
p_r{2} = R_y.*ear_r_y + R;
p_r{3} = R_z.*ear_r_z;

%% Left Ear
[ear_l_x, ear_l_y, ear_l_z] = sph2cart(Azi, Ele, 1);
p_l{1} = R_x.*ear_l_x;
p_l{2} = R_y.*ear_l_y - R;
p_l{3} = R_z.*ear_l_z;

end