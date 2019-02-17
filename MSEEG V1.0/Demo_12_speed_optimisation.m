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

%% Initialisation
clear; close all; clc;
addpath(genpath('Toolbox'));
addpath(genpath('Data\Optimisation'));
addpath(genpath('Data\Power Maps'));
addpath(genpath('Data\EEG'));

%% Parameters
M    = 7;        % Total number of patients
fs   = 32;       % Sampling frequency
T    = 256;      % Segment length in samples

%% Loading Database information
load('Ratios');
load('Data\Optimisation\Speed\Optimal_Speeds.mat');

%% Main
C = zeros(M,1);
E = zeros(M,1);
for m = 1:7
    % Loading Data
    load(['Power Maps\Patient_' num2str(m) '_power.mat'],'P_eeg');
    load(['Data\Optimisation\Speed\Patient_' num2str(m) '.mat'],'X_best');
    % Computing simulated power maps
    P_sim = zeros(size(P_eeg));
    for nseg = 1:size(P_eeg,2)
        st = 1 + (nseg-1)*T;
        fi = nseg*T;
        P_sim(:,nseg) = sum(X_best(:,st:fi).^2,2)./T;
    end
    % Computing error and correlation
    C(m,1) = 100*abs(corr(P_sim(:),P_eeg(:)));
    E(m,1) = 100*norm(P_eeg - P_sim)./norm(P_eeg);
    fprintf('\t  %1.0f\t\t\t%6.4f\t\t%6.4f\t\t%6.4f\n',m,V(m),C(m),E(m));
end
sqrt(sum((C-mean(sum(R(:,1).*C))).^2)/6)

fprintf('\twAVG\t\t\t\t\t%6.4f\t\t%6.4f\n',sum(R(:,1).*C),sum(R(:,1).*E));
fprintf('\tSTD\t\t\t\t\t\t%6.4f\t\t%6.4f\n',sqrt(sum((C-mean(sum(R(:,1).*C))).^2)/6),sqrt(sum((E-mean(sum(R(:,1).*E))).^2)/6));

