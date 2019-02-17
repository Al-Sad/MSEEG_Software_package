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
Q = 100;      % Sample size of noise level
T = 256;      % Number of samples in one segment
rng(1);       % Ensuring reproducibility

%% Loading noise estimates
load('Noise_Estimation')

%% Main
N_opt_s = [];
N_opt_b = [];
NL    = linspace(1,100,Q);
for m = 1:M
    % Loading Data
    disp(['Patient : ' num2str(m)])
    load(['Patient_' num2str(m)],'X_eeg','X_sim','mask');
    
    % Extracting background and seizure epochs
    mask = logical(mask);
    X_sim_s  = X_sim(:,mask);
    X_real_s = X_eeg(:,mask);
    X_sim_b  = X_sim(:,~mask);
    X_real_b = X_eeg(:,~mask);
    
    % Initialisation
    Y_sim_s = zeros(size(X_sim_s));
    Y_sim_b = zeros(size(X_sim_b));
    Ls      = length(X_real_s)/T;
    Lb      = length(X_real_b)/T;
    Er_s    = zeros(Ls,Q);
    Er_b    = zeros(Lb,Q);

    % Optimising seizures NL
    disp('Optimising seizure noise level');
    Is = zeros(1,Ls);
    for i = 1:Ls
        % Segmenting
        st = 1 + (i-1)*T;
        fi = i*T;
        X_real_seg_s = X_real_s(:,st:fi);
        X_sim_seg_s  = X_sim_s(:,st:fi);
        % Optimising randomness
        Cr = abs(corr(X_real_seg_s',X_real_seg_s'));
        for q = 1:Q
            Ys_seg = power_div(X_sim_seg_s,NL(q),D_s,m_s,s_s,n_s);
            Cy = abs(corr(Ys_seg',Ys_seg'));
            Er_s(i,q) = norm(Cr - Cy)/norm(Cr);
        end
        Is(1,i) = find(Er_s(i,:) == min(Er_s(i,:)));
        disp(100*i/Ls)
    end
    N_opt_s = [N_opt_s NL(Is)];
    
    % Optimising background NL
    disp('Optimising background noise level');
    Ib = zeros(1,Lb);
    for i = 1:Lb
        % Segmenting
        st = 1 + (i-1)*T;
        fi = i*T;
        X_real_seg_b = X_real_b(:,st:fi);
        X_sim_seg_b  = X_sim_b(:,st:fi);
        % Optimising randomness
        Cr = abs(corr(X_real_seg_b',X_real_seg_b'));
        for q = 1:Q
            Yb_seg = power_div(X_sim_seg_b,NL(q),D_b,m_b,s_b,n_b);
            Cy = abs(corr(Yb_seg',Yb_seg'));
            Er_b(i,q) = norm(Cr - Cy)/norm(Cr);
        end
        Ib(1,i) = find(Er_b(i,:) == min(Er_b(i,:)));
        disp(100*i/Lb)
    end
    N_opt_b = [N_opt_b NL(Ib)];
end

%% Saving
save('Data\Optimisation\Noise\Noise_Level','N_opt_s','N_opt_b');
