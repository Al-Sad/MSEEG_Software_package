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
rng(1);          % Ensuring reproducibility

%% Loading optimal intensities and speeds
load('Intensity\Optimal_Intensities');
load('Speed\Optimal_Speeds');

%% Main
for m = 1:M
    % Loading Data
    disp(['Patient # ' num2str(m)]);
    load(['Power Maps\Patient_' num2str(m) '_power.mat'],'P_eeg','mask','mask_seg');
    load(['EEG\Patient_' num2str(m) '_Data_corrected.mat'],'clean_EEG','Channel_label');
    load(['Positions\Optimal Positions\Patient_' num2str(m) '.mat']);
    % Reshaping Real EEG
    X_eeg_seg = clean_EEG;
    [Nseg, ch_n, N] = size(X_eeg_seg);
    X_eeg = zeros(ch_n, N*Nseg);
    for nseg = 1:Nseg
        st = 1 + (nseg-1)*T;
        fi = nseg*T;
        X_eeg(:,st:fi) = X_eeg_seg(nseg,:,:);
    end
    % Generating simulated multi-sensor EEG
    X_sim = multi_sensor_EEG(mask, fs, opt_pos_back, opt_pos_seiz, mP_b(m), mP_s(m), V(m));
    % Computing simulated power maps
    P_sim = zeros(size(P_eeg));
    X_sim_seg = zeros(size(X_eeg_seg));
    for nseg = 1:Nseg
        st = 1 + (nseg-1)*T;
        fi = nseg*T;
        P_sim(:,nseg) = sum(X_sim(:,st:fi).^2,2)./T;
        X_sim_seg(nseg,:,:) = X_sim(:,st:fi);
    end
    % Saving
    save(['Data\Simulated EEG\Patient_' num2str(m)],'X_eeg','X_sim','X_sim_seg','P_eeg','P_sim','fs','Channel_label','mask','mask_seg');
end
