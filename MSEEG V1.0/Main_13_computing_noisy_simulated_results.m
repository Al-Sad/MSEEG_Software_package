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
M = 7;        % Total number of patients
T = 256;      % Number of samples in one segment
rng(1);       % Ensuring reproducibility

%% Loading noise signal
load('Noise_Estimation')
DD_s = D_s;
mm_s = m_s;
ss_s = s_s;
nn_s = n_s;
DD_b = D_b;
mm_b = m_b;
ss_b = s_b;
nn_b = n_b;
load('Noise_Level_distribution');

%% Main
for m = 1:M
    % Loading Data
    disp(['Patient : ' num2str(m)])
    load(['Patient_' num2str(m)],'X_sim','X_sim_seg','mask','mask_seg');
    
    % Extracting segments
    mask     = logical(mask);
    mask_seg = logical(mask_seg);
    Xs = X_sim(:,mask);
    Xb = X_sim(:,~mask);
    [ch_n, N] = size(X_sim);
    Ns = length(Xs);
    Nb = length(Xb);
    
    % Initialisation
    X_sim_noisy     = X_sim;
    X_sim_seg_noisy = X_sim_seg;
    tmp_s = zeros(Ns/T, ch_n, T);
    tmp_b = zeros(Nb/T, ch_n, T);
    Xs_r  = zeros(size(Xs));
    Xb_r  = zeros(size(Xb));
    
    % Segmenting seizures and adding noise
    for i = 1:Ns/T
        st = 1 + (i-1)*T;
        fi = i*T;
        NL = abs(random(D_s,M_s,S_s));
        if(NL >= 100), NL = 100; end
        tmp_s(i,:,:) = real(power_div(Xs(:,st:fi),NL,DD_s,mm_s,ss_s,nn_s));
        Xs_r(:,st:fi) = tmp_s(i,:,:);
    end
    % Segmenting background and adding noise
    for i = 1:Nb/T
        st = 1 + (i-1)*T;
        fi = i*T;
        NL = abs(random(D_b,M_b,S_b));
        if(NL >= 100), NL = 100; end
        tmp_b(i,:,:) = real(power_div(Xb(:,st:fi),NL,DD_b,mm_b,ss_b,nn_b));
        Xb_r(:,st:fi) = tmp_b(i,:,:);
    end
    
    % Appending noisy signals
    X_sim_noisy(:,mask)            = Xs_r;
    X_sim_seg_noisy(mask_seg,:,:)  = tmp_s;
    X_sim_noisy(:,~mask)           = Xb_r;
    X_sim_seg_noisy(~mask_seg,:,:) = tmp_b;

    % Saving
    save(['Data\Simulated EEG\Patient_' num2str(m)],'X_sim_noisy','X_sim_seg_noisy','-append');
end