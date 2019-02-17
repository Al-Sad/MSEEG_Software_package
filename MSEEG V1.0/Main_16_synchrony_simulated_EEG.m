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
M    = 7;     % Total number of patients
Nfft = 256;   % Total number of frequency bins

%% Main
for m = 1:M
    disp(['Patient No. ' num2str(m)])
    % Loading Real EEG
    load(['Patient_' num2str(m)],'X_sim_seg_noisy','mask_seg','Channel_label');
    Ne = size(X_sim_seg_noisy,1);
    Nb   = size(X_sim_seg_noisy,1)-sum(mask_seg);
    Ns   = sum(mask_seg);
    ch_n = size(X_sim_seg_noisy,2);
    N    = size(X_sim_seg_noisy,3);
    EEG_S = zeros(Ns,ch_n,N);
    EEG_B = zeros(Nb,ch_n,N);
    cnt1 = 0; cnt2 = 0;
    for ne = 1:Ne
        if(mask_seg(ne))
            cnt1 = cnt1 + 1;
            EEG_S(cnt1,:,:) = X_sim_seg_noisy(ne,:,:);
        else
            cnt2 = cnt2 + 1;
            EEG_B(cnt2,:,:) = X_sim_seg_noisy(ne,:,:);
        end
    end
    CB  = cell(1,Nb);
    CS  = cell(1,Ns);
    tmp = [];
    % Causality analysis of background EEG
    disp('Causality analysis of background EEG ...');
    for nb = 1:Nb
        tmp(:,:) = EEG_B(nb,:,:);
        D = mtfd(tmp, 'ckd', 1, 0.02, 0.04, Nfft);
        for i = 1:ch_n
            D{i,i} = real(D{i,i});
        end
        CB{1,nb} = PLV(D, 'extended'); clear D;
        disp(100*nb/Nb);
    end
    tmp = [];
    % Causality analysis of seizure EEG
    disp('Causality analysis of seizure EEG ...');
    for ns = 1:Ns
        tmp(:,:) = EEG_S(ns,:,:);
        D = mtfd(tmp, 'ckd', 1, 0.02, 0.04, Nfft);
        for i = 1:ch_n
            D{i,i} = real(D{i,i});
        end
        CS{1,ns} = PLV(D, 'extended'); clear D;
        disp(100*ns/Ns);
    end
    % Saving
    disp('Saving ...')
    save(['Data\Synchrony\simulated_EEG_' num2str(m) '_synch.mat'],'CB','CS','Channel_label');
end
