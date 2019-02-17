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
%                  Sphere cross-section polygon model
%
%  Syntax  : p = sphere_cross(R, azi_s, azi_f, N, M)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R        : Sphere radius in cm.
% azi_s    : Start of azimuth spanning angle.
% azi_f    : End of azimuth spanning angle.
% N        : Meshing resolution of the elevation angle.
% M        : Meshing resolution of the azimuth angle.
%
% <OUTPUTs>
% p        : Sphere cross-section polygon in cartesian coordinates.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% for i = 1:0.01:3
%     p = sphere_cross(i, -pi/2, pi, 32, 32);
%     surf(p{1}, p{2}, p{3},'Linestyle','-','EdgeAlpha',0.2); hold on;
% end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function p = sphere_cross(R, azi_s, azi_f, N, M)
azi = linspace(azi_s, azi_f, M);
ele = linspace(-pi/2, pi/2, N);
[Azi, Ele] = meshgrid(azi, ele);
[p{1}, p{2}, p{3}] = sph2cart(Azi, Ele, R);
end