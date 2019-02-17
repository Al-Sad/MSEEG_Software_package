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
%                         Noise level computation
%
% Syntax   : [y, Py, Px, Pn] = power_div(x,NL,varargin)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% x        : Clean EEG segment.
% NL       : Noise level in percentage (0 to 100).
% varargin : This is a dynamic input argument that has to be supplied in
%            the following sequence:
%            varargin{1} : PDF with three parameters name, see "help random".
%            varargin{2} : PDF first parameter.
%            varargin{3} : PDF second parameter.
%            varargin{4} : PDF third parameter.
%
% <OUTPUTs>
% y        : Noisy EEG segment.
% Py       : Power of generated noisy EEG.
% Px       : Power of input clean EEG.
% Pn       : Power of generated random noise.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% x = EEG_seiz(256,20);
% [y, Py, Px, Pn] = power_div(x,40,'tLocationScale',0,13,6);
% figure; plot(x); hold on; plot(y,'r'); legend('clean','noisy');
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [y, Py, Px, Pn] = power_div(x,NL,varargin)
%% Initialisation
[M, N] = size(x);
y      = zeros(size(x));
Py     = zeros(M,1);
Px     = zeros(M,1);
Pn     = zeros(M,1);

%% Main
for m = 1:M
    px  = sum(x(m,:).^2)./N;
    x(m,:) = sqrt(1-NL/100).*x(m,:);
    n = random(varargin{1},varargin{2},varargin{3},varargin{4},1,N);
    pn  = sum(n.^2)./N;
    n   = n ./ sqrt(pn);
    tmp = n;
    
    %% Secant Method for Beta finding
    xn_2 = px;
    xn_1 = (1-NL/100)*px;
    f2 = px - power_fun(x(m,:), tmp, xn_2);
    f1 = px - power_fun(x(m,:), tmp, xn_1);
    cost = Inf;
    iter = 0;
    while(cost>1e-6)
        iter = iter + 1;
        x_n = abs(xn_1 - f1*((xn_1 - xn_2)/(f1 - f2)));
        xn_2 = xn_1;
        xn_1 = x_n;
        f2 = px - power_fun(x(m,:), tmp, xn_2);
        f1 = px - power_fun(x(m,:), tmp, xn_1);
        cost = abs(px - power_fun(x(m,:), tmp, x_n));
    end
    
    %% Adding noise
    Px(m,:)  = sum(x(m,:).^2)./N;
    Pn(m,:)  = x_n;
    [Py(m,1), y(m,:)] = power_fun(x(m,:), tmp, x_n);
end
end