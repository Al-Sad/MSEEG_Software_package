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

%% Parameters
M    = 7;        % Total number of patients
fs   = 32;       % Sampling frequency
T    = 256;      % Segment length in samples
Q    = 1000;     % Number of iterations per patient

%% Loading optimal intensities and speeds
load('Data\Optimisation\Intensity\Optimal_Intensities');

%% Main 1
CC = zeros(M,1);
EE = zeros(M,1);
V  = zeros(M,1);
fprintf('\t  #\t\t\t  v\t\t\t Corr\t\t Error\n');
for m = 1:M
    load(['Data\Optimisation\Speed\Patient_' num2str(m)]);
    P = 0.5.*(100 + C - E);
    PP = max(max(P));
    [Ir, Ic] = find(P==PP);
    CC(m) = C(Ir,Ic); EE(m) = E(Ir,Ic); V(m) = v(1,Ic);
end
save('Data\Optimisation\Speed\Optimal_Speeds','V');

%% Main 2
for m = 1:M
    rng(1);          % Ensure reproducibility
    c = 0;
    % Loading Data
    disp(['Patient # ' num2str(m)]);
    load(['Data\Optimisation\Positions\Optimal Positions\Patient_' num2str(m)]);
    load(['Power Maps\Patient_' num2str(m) '_power.mat'],'P_eeg','mask');
    E = zeros(Q,1);
    C = zeros(Q,1);
    P = zeros(Q,1);
    for q = 1:Q
        % Generating simulated multi-sensor EEG
        X = multi_sensor_EEG(mask, fs, opt_pos_back, opt_pos_seiz, mP_b(m), mP_s(m), V(m));
        % Computing simulated power maps
        P_sim = zeros(size(P_eeg));
        for nseg = 1:size(P_eeg,2)
            st = 1 + (nseg-1)*T;
            fi = nseg*T;
            P_sim(:,nseg) = sum(X(:,st:fi).^2,2)./T;
        end
        % Computing error and correlation
        E(q,1) = 100*norm(P_eeg - P_sim)./norm(P_eeg);
        C(q,1) = 100*corr(P_sim(:),P_eeg(:));
        P(q,1) = 0.5.*(100 + C(q,1) - E(q,1));
        
        % Save best matching signal
        if(P(q,1)>c)
            c = P(q,1);
            X_best = X;
        end
        disp(100*q/Q)
    end
    % Saving
    save(['Data\Optimisation\Speed\Patient_' num2str(m)],'X_best','-append');
end
