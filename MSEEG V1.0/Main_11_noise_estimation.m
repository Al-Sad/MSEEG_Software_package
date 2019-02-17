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
addpath(genpath('Data\Simulated EEG'));
addpath(genpath('Data\Optimisation\Noise'));

%% Parameters
M   = 7;    % Total number of patients
bin = 256;  % Total number of bins in a histogram

%% Main
G_s = [];
G_b = [];
for m = 1:M
    % Loading Data
    disp(['Patient # ' num2str(m)]);
    load(['Patient_' num2str(m)],'X_eeg','mask');
    Xs_real = X_eeg(:,logical(mask));
    Xb_real = X_eeg(:,~logical(mask));
    G_s = [G_s Xs_real(:)'];
    G_b = [G_b Xb_real(:)'];
end

%% PDF Estimation
[H_s, x_s] = hist(G_s,bin);
H_s = H_s ./ trapz(x_s,H_s);
[H_b, x_b] = hist(G_b,bin);
H_b = H_b ./ trapz(x_b,H_b);
L_s = fitdist(G_s','tLocationScale');
L_b = fitdist(G_b','tLocationScale');
x = -70:0.1:70;
P_s = pdf('tLocationScale',x,L_s.mu,L_s.sigma,L_s.nu);
P_b = pdf('tLocationScale',x,L_b.mu,L_b.sigma,L_b.nu);

%% Saving
D_s = 'tlocationscale';
m_s = L_s.mu;
s_s = L_s.sigma;
n_s = L_s.nu;
D_b = 'tlocationscale';
m_b = L_b.mu;
s_b = L_b.sigma;
n_b = L_b.nu;
save('Data\Optimisation\Noise\Noise_Estimation','H_s','H_b','x_s','x_b',...
    'D_s','m_s','s_s','n_s','D_b','m_b','s_b','n_b');