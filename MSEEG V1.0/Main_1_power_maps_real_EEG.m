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
addpath(genpath('Data\EEG'));

%% Parameters
M   = 7;     % Total number of patients
thr = 0.998; % PDF threshold of EEG power

%% Calculating Power Threshold
disp('Calculating Power Threshold')
P = [];
for m = 1:M
    % Loading Real EEG
    load(['Patient_' num2str(m) '_Data.mat'],'clean_EEG');
    % Multi-Channel Power
    [Nseg, ch_n, N] = size(clean_EEG);
    Peeg = zeros(ch_n, Nseg);
    for i = 1:Nseg
        tmp(:,:)   = clean_EEG(i,:,:);
        Peeg(:,i) = sum(tmp.^2,2)./N;
    end
    P = [P Peeg(:)'];
end
% Defining 99.85% CI of Powers
[tmp, x] = hist(P,length(P));
for i = 1:length(tmp)
    PD = sum(tmp(1:i))/length(P);
    if(PD >= thr)
        Cutoff = x(i);
        break
    end
end

%% Computing and Correcting Power Maps
disp('Computing and Correcting Power Maps')
for m = 1:M
    tmp = [];
    disp(['Patient No. ' num2str(m)])
    % Loading Real EEG
    load(['Patient_' num2str(m) '_Data.mat'],'clean_EEG','EEG_S','EEG_B','EEG_mask','Channel_label','fs');
    % Computing Power maps
    [Nseg, ch_n, N] = size(clean_EEG);
    Peeg = zeros(ch_n, Nseg);
    for i = 1:Nseg
        tmp(:,:)   = clean_EEG(i,:,:);
        Peeg(:,i) = sum(tmp.^2,2)./N;
    end
    [Ns, ~, N] = size(EEG_S);
    Ps = zeros(ch_n, Ns);
    for i = 1:Ns
        tmp(:,:) = EEG_S(i,:,:);
        Ps(:,i) = sum(tmp.^2,2)./N;
    end
    [Nb, ~, N] = size(EEG_B);
    Pb = zeros(ch_n, Nb);
    for i = 1:Nb
        tmp(:,:) = EEG_B(i,:,:);
        Pb(:,i) = sum(tmp.^2,2)./N;
    end
    % Correcting Power maps
    P_s   = Ps.*(Ps <= Cutoff) + Cutoff.*(Ps > Cutoff);
    P_b   = Pb.*(Pb <= Cutoff) + Cutoff.*(Pb > Cutoff);
    P_eeg = Peeg.*(Peeg <= Cutoff) + Cutoff.*(Peeg > Cutoff);
    % Correcting EEG data
    Ac_s = sqrt(P_s./Ps);
    Ac_b = sqrt(P_b./Pb);
    Ac_e = sqrt(P_eeg./Peeg);
    for i = 1:Ns
        tmp = []; tmp(:,:) = EEG_S(i,:,:);
        for j = 1:ch_n
            EEG_S(i,j,:) = Ac_s(j,i).*tmp(j,:);
        end
    end
    for i = 1:Nb
        tmp = []; tmp(:,:) = EEG_B(i,:,:);
        for j = 1:ch_n
            EEG_B(i,j,:) = Ac_b(j,i).*tmp(j,:);
        end
    end
    for i = 1:Nseg
        tmp = []; tmp(:,:) = clean_EEG(i,:,:);
        for j = 1:ch_n
            clean_EEG(i,j,:) = Ac_e(j,i).*tmp(j,:);
        end
    end
    % Saving
    disp('Saving ...')
    save(['Data\EEG\Patient_' num2str(m) '_Data_corrected.mat'],...
        'clean_EEG','EEG_S','EEG_B','EEG_mask','Channel_label','fs');
    save(['Data\Power Maps\Patient_' num2str(m) '_power.mat'],...
        'P_eeg','P_s','P_b','Channel_label');
end

%% Rehaping and saving EEG masks
for m = 1:M
    mask = [];
    load(['Patient_' num2str(m) '_Data_corrected.mat'],'EEG_mask');
    Ne = length(EEG_mask);
    mask_seg = EEG_mask;
    for ne = 1:Ne
        mask = [mask repmat((EEG_mask(ne)==1),[1 256])];
    end
    save(['Data\Power Maps\Patient_' num2str(m) '_power.mat'],'mask','mask_seg','-append');
end
