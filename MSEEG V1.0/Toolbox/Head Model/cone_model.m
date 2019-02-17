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
%                          Nose polygon model
%
%  Syntax  : p = cone_model(R, azi_s, azi_f, ele, L, N, M)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R        : Head radius in cm to center the cone.
% azi_s    : Start of azimuth spanning angle.
% azi_f    : End of azimuth spanning angle.
% ele      : Elevation angle for orientation.
% L        : Cone height in cm.
% N        : Meshing resolution along the cone height.
% M        : Meshing resolution along the cone base plane.
%
% <OUTPUTs>
% p        : Cone polygon in cartesian coordinates.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% p = cone_model(5, 0, 2*pi, pi/3, 2, 32, 32);
% surf(p{1},p{2},p{3});
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function p = cone_model(R, azi_s, azi_f, ele, L, N, M)
azi = linspace(azi_s, azi_f, M);
r   = linspace(0, L, N);
[Azi, RR] = meshgrid(azi,r);
p{1} = RR.*sin(ele) - R - L*sin(ele) + 0.1;
p{2} = RR.*cos(ele).*sin(Azi);
p{3} = RR.*cos(ele).*cos(Azi);
end