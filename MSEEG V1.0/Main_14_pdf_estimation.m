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

%% Parameters
M   = 7;    % Total number of patients
bin = 256;  % Total number of bins in a histogram

%% Main
G_sim_s = [];
G_nos_s = [];
G_sim_b = [];
G_nos_b = [];
for m = 1:M
    % Loading Data
    disp(['Patient # ' num2str(m)]);
    load(['Patient_' num2str(m)],'X_sim','X_sim_noisy','mask');
    Xs_s = X_sim(:,logical(mask));
    Xn_s = X_sim_noisy(:,logical(mask));
    Xs_b = X_sim(:,~logical(mask));
    Xn_b = X_sim_noisy(:,~logical(mask));
    G_sim_s = [G_sim_s Xs_s(:)']; clear  Xs_s;
    G_nos_s = [G_nos_s Xn_s(:)']; clear  Xn_s;
    G_sim_b = [G_sim_b Xs_b(:)']; clear  Xs_b;
    G_nos_b = [G_nos_b Xn_b(:)']; clear  Xn_b;
end

%% PDF Estimation
[H_sim_s, x_sim_s] = hist(G_sim_s,bin);
H_sim_s = H_sim_s ./ trapz(x_sim_s,H_sim_s);
[H_nos_s, x_nos_s] = hist(G_nos_s,bin);
H_nos_s = H_nos_s ./ trapz(x_nos_s,H_nos_s);
[H_sim_b, x_sim_b] = hist(G_sim_b,bin);
H_sim_b = H_sim_b ./ trapz(x_sim_b,H_sim_b);
[H_nos_b, x_nos_b] = hist(G_nos_b,bin);
H_nos_b = H_nos_b ./ trapz(x_nos_b,H_nos_b);
L_sim_s = fitdist(G_sim_s','tLocationScale');
L_nos_s = fitdist(G_nos_s','tLocationScale');
L_sim_b = fitdist(G_sim_b','tLocationScale');
L_nos_b = fitdist(G_nos_b','tLocationScale');

%% Saving
D_sim_s = 'tlocationscale';
m_sim_s = L_sim_s.mu;
s_sim_s = L_sim_s.sigma;
n_sim_s = L_sim_s.nu;
D_nos_s = 'tlocationscale';
m_nos_s = L_nos_s.mu;
s_nos_s = L_nos_s.sigma;
n_nos_s = L_nos_s.nu;
D_sim_b = 'tlocationscale';
m_sim_b = L_sim_b.mu;
s_sim_b = L_sim_b.sigma;
n_sim_b = L_sim_b.nu;
D_nos_b = 'tlocationscale';
m_nos_b = L_nos_b.mu;
s_nos_b = L_nos_b.sigma;
n_nos_b = L_nos_b.nu;
save('Data\Optimisation\Noise\Clean_PDF','H_sim_s','H_sim_b','x_sim_s','x_sim_b',...
    'D_sim_s','m_sim_s','s_sim_s','n_sim_s','D_sim_b','m_sim_b','s_sim_b','n_sim_b');
save('Data\Optimisation\Noise\Noisy_PDF','H_nos_s','H_nos_b','x_nos_s','x_nos_b',...
    'D_nos_s','m_nos_s','s_nos_s','n_nos_s','D_nos_b','m_nos_b','s_nos_b','n_nos_b');