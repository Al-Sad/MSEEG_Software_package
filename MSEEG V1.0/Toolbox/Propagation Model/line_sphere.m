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
%                   Line-and-sphere intersection method
%
% Syntax  : pl = line_sphere(event_p, elec_p, R)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% event_p : Single EEG event location within the brain sphere.
% elec_p  : Electrodes location on the scalp sphere.
% R       : Sphere radius in cm.
%
% <OUTPUTs>
% pl      : Cartesian positions of the projected electrodes on the sphere
%           with radius R.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% event_p = event_pos(1);
% elec_p = electrode_10_20(5.95);
% pl = line_sphere(event_p, elec_p, 5);
% disp(pl);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function pl = line_sphere(event_p, elec_p, R)
%% Parameters
tol  = 1e-20;
[ch_n,~] = size(elec_p);        % Number of electrodes

%% Electrode-Event Relative Distances
dx = elec_p(:,1) - repmat(event_p(1), ch_n, 1);
dy = elec_p(:,2) - repmat(event_p(2), ch_n, 1);
dz = elec_p(:,3) - repmat(event_p(3), ch_n, 1);
dd = [dx dy dz];
[Azm, Ele, ~]   = cart2sph(dx, dy, dz);
pl = zeros(ch_n,3);

for i = 1:ch_n
    %% Computing Sphere equation
    % coefficients
    a = sum(dd(i,:) .* dd(i,:), 2);
    b = 2*sum(event_p.*dd(i,:), 2);
    c = sum(event_p.*event_p, 2) - R^2;
    % solve equation
    delta = b.*b - 4*a.*c;
    
    %% Finding Intersections
    if delta > tol
        % delta positive: find two roots of second order equation
        u1 = (-b -sqrt(delta)) / 2 / a;
        u2 = (-b +sqrt(delta)) / 2 / a;
        % convert into 3D coordinate
        point = [event_p+u1*dd(i,:) ; event_p+u2*dd(i,:)];
        
    elseif abs(delta) < tol
        % delta around zero: find unique root, and convert to 3D coord.
        u = -b/2./a;
        point = repmat(event_p + u*dd(i,:),2,1);
    else
        % delta negative: no solution
        point = ones(2, 3);
        point(:) = NaN;
    end
    
    %% Ommiting irrelevant intersecction points
    [Azm1, Ele1, ~] = cart2sph(point(1,1)-event_p(1,1), point(1,2)-event_p(1,2), point(1,3)-event_p(1,3));
    [Azm2, Ele2, ~] = cart2sph(point(2,1)-event_p(1,1), point(2,2)-event_p(1,2), point(2,3)-event_p(1,3));
    
    if( (abs(Azm1-Azm(i))<1e-6) && (abs(Ele1-Ele(i))<1e-6) )
        pl(i,:) = point(1,:);
    elseif( (abs(Azm2-Azm(i))<1e-6) && (abs(Ele2-Ele(i))<1e-6) )
        pl(i,:) = point(2,:);
    end
end
end