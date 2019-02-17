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
addpath(genpath('Data\Optimisation\Intensity'));
addpath(genpath('Data\EEG'));

%% Parameters
M = 7; % Total number of patients

%% Loading Database information
load('Ratios');

%% Initialise
mP_b = zeros(M,1);
mP_s = zeros(M,1);
CC_b = zeros(M,1);
CC_s = zeros(M,1);
EE_b = zeros(M,1);
EE_s = zeros(M,1);

%% Display results
fprintf('\t\t\t\t\t\t Background\n');
fprintf('\t  #\t\t\t  Io\t\t Corr\t\t Error\n');
for m = 1:M
    load(['Patient_' num2str(m) '_background.mat']);
    PP = 0.5.*(100 + Corr - Error);
    [~, I] = max(PP);
    mP_b(m) = mP(I);
    CC_b(m) = Corr(I);
    EE_b(m) = Error(I);
    fprintf('\t  %1.0f\t\t\t%6.1f\t\t%6.4f\t\t%6.4f\n',m,mP_b(m),CC_b(m),EE_b(m));
end
fprintf('\twAVG\t\t\t\t\t%6.4f\t\t%6.4f\n',sum(R(:,3).*CC_b),sum(R(:,3).*EE_b));

fprintf('\n\n\t\t\t\t\t\t Seizure\n');
fprintf('\t  #\t\t\t  Io\t\t Corr\t\t Error\n');
for m = 1:M
    load(['Patient_' num2str(m) '_seizure.mat']);
    PP = 0.5.*(100 + Corr - Error);
    [~, I] = max(PP);
    mP_s(m) = mP(I);
    CC_s(m) = Corr(I);
    EE_s(m) = Error(I);
    fprintf('\t  %1.0f\t\t\t%6.1f\t\t%6.4f\t\t%6.4f\n',m,mP_s(m),CC_s(m),EE_s(m));
end
fprintf('\twAVG\t\t\t\t\t%6.4f\t\t%6.4f\n',sum(R(:,2).*CC_s),sum(R(:,2).*EE_s));

%% Save Results
save('Data\Optimisation\Intensity\Optimal_Intensities','mP_b','mP_s');