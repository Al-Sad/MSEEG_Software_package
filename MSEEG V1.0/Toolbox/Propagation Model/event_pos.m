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
%                   Random EEG source position generation
%
% Syntax : p = event_pos(N)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N      : Number of source signals inside the brain sphere.
%
% <OUTPUTs>
% p      : Random source signals position in 3D cartesian coordinates.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% p = event_pos(5);
% disp(p);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function p = event_pos(N)
%% Parameters
Scalp_R  = 5.95;        % Head Radius in cm
th1      = 0.29;        % Scalp thickness in cm
th2      = 0.6;         % Skull thickness in cm
th3      = 0.3;         % CSF thickness in cm
brain_R  = Scalp_R - th1 - th2 - th3;

%% Main
a = 0; b = brain_R-0.01;
R = (rand(N,1)*(b^3-a^3)+a^3).^(1/3);
phi1 = acos(rand(N,1));
th1  = 2*pi*rand(N,1);
x    = R.*sin(phi1).*sin(th1);
y    = R.*sin(phi1).*cos(th1);
z    = R.*cos(phi1);
p    = [x, y, z];
end