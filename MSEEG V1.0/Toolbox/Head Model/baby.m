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
%                     Baby head polygon normlaization
%
%  Syntax  : [Faces,Vertices] = baby(R)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% R        : Head radius in cm.
%
% <OUTPUTs>
% Faces    : Polygon face values.
% Vertices : Polygon vertex values.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% [F,V] = baby(4);
% trisurf(F,V(:,1),V(:,2),V(:,3));
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [Faces, Vertices] = baby(R)
%% Loading Head polygon
load('head.mat','Vertices','Faces');
%% Normalization according to the specified head radius R
Vertices(:,1) = R.*(-1).*Vertices(:,1)./max(abs(Vertices(:,1)));
Vertices(:,2) = R.*(-1).*Vertices(:,2)./max(abs(Vertices(:,2)));
Vertices(:,3) = R.*Vertices(:,3)./max(abs(Vertices(:,3)));
end
