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
%                       Adding multi-sensor noise
%
% Syntax : [py, y] = power_fun(x, tmp, alpha)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% x      : Clean EEG segment.
% tmp    : Normalized random noise.
% alpha  : Noise power factor (0 to 1).
%
% <OUTPUTs>
% py     : Power of generated noisy EEG.
% y      : Noisy EEG segment.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% x = EEG_seiz(256,20);
% n = random('tLocationScale',0,13,6,[1 256]);
% n = n./sqrt(sum(n.^2)./256);
% [py, y] = power_fun(x, n, 0.3);
% figure; plot(x); hold on; plot(y,'r'); legend('clean','noisy');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [py, y] = power_fun(x, tmp, alpha)
n = sqrt(alpha).*tmp;
y  = x + n;
py = sum(y.^2)./length(x);
end