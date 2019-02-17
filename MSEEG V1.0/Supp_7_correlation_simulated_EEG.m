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
addpath(genpath('Data\Correlation'));
addpath(genpath('Data\EEG'));

%% Parameters
M    = 7;     % Total number of patients
Nfft = 256;   % Total number of frequency bins

%% Main
load('Ratios');
out_t  = zeros(M,6);
out_tf = zeros(M,6);
for m = 1:M
    disp(['Patient No. ' num2str(m)]);
    % Loading Real EEG
    load(['simulated_EEG_' num2str(m) '_corr.mat'],'CB_t','CS_t','CB_tf','CS_tf');
    ch_n = length(CB_t{1,1});
    % Averaging
    disp('Averaging all patients ...');
    Nb = length(CB_t);
    Corr_B_t  = zeros(1,Nb*(ch_n^2-ch_n)/2);
    Corr_B_tf = zeros(1,Nb*(ch_n^2-ch_n)/2);
    for nb = 1:Nb
        st = 1 + (nb-1)*(ch_n^2-ch_n)/2;
        fi = nb*(ch_n^2-ch_n)/2;
        Corr_B_t(1,st:fi)  = nonzeros(triu(CB_t{1,nb},1));
        Corr_B_tf(1,st:fi) = nonzeros(triu(CB_tf{1,nb},1));
    end
    Ns = length(CS_t);
    Corr_S_t  = zeros(1,Ns*(ch_n^2-ch_n)/2);
    Corr_S_tf = zeros(1,Ns*(ch_n^2-ch_n)/2);
    for ns = 1:Ns
        st = 1 + (ns-1)*(ch_n^2-ch_n)/2;
        fi = ns*(ch_n^2-ch_n)/2;
        Corr_S_t(1,st:fi)  = nonzeros(triu(CS_t{1,ns},1));
        Corr_S_tf(1,st:fi) = nonzeros(triu(CS_tf{1,ns},1));
    end
    % Arranging output
    out_t(m,1:3)  = [min(Corr_B_t) mean(Corr_B_t) max(Corr_B_t)];
    out_t(m,4:6)  = [min(Corr_S_t) mean(Corr_S_t) max(Corr_S_t)];
    out_tf(m,1:3) = [min(Corr_B_tf) mean(Corr_B_tf) max(Corr_B_tf)];
    out_tf(m,4:6) = [min(Corr_S_tf) mean(Corr_S_tf) max(Corr_S_tf)];
end
M_t  = zeros(1,6);
M_tf = zeros(1,6);
for m = 1:M
    M_t  = [R(m,3).*out_t(m,1:3) R(m,2).*out_t(m,4:6)] + M_t;
    M_tf = [R(m,3).*out_tf(m,1:3) R(m,2).*out_tf(m,4:6)] + M_tf;
end

%% Dispaly Results
disp('Time domain correlation results');
disp(out_t); disp(M_t);
disp('Time-frequency correlation results');
disp(out_tf); disp(M_tf);
