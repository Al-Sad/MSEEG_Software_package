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
addpath(genpath('Data\Synchrony'));
addpath(genpath('Data\EEG'));

%% Parameters
M = 7;     % Total number of patients

%% Main
load('Ratios');
out  = zeros(M,6);
for m = 1:M
    disp(['Patient No. ' num2str(m)]);
    % Loading Real EEG
    load(['Patient_' num2str(m) '_synch.mat'],'CB','CS');
    ch_n = length(CB{1,1});
    % Averaging
    disp('Averaging all results ...');
    Nb = length(CB);
    Synch_B  = zeros(1,Nb*(ch_n^2-ch_n)/2);
    for nb = 1:Nb
        st = 1 + (nb-1)*(ch_n^2-ch_n)/2;
        fi = nb*(ch_n^2-ch_n)/2;
        Synch_B(1,st:fi) = nonzeros(triu(abs(CB{1,nb}),1));
    end
    Ns = length(CS);
    Synch_S  = zeros(1,Ns*(ch_n^2-ch_n)/2);
    for ns = 1:Ns
        st = 1 + (ns-1)*(ch_n^2-ch_n)/2;
        fi = ns*(ch_n^2-ch_n)/2;
        Synch_S(1,st:fi) = nonzeros(triu(abs(CS{1,ns}),1));
    end
    % Arranging output
    out(m,1:3) = [min(Synch_B) mean(Synch_B) max(Synch_B)];
    out(m,4:6) = [min(Synch_S) mean(Synch_S) max(Synch_S)];
end
M_tf = zeros(1,6);
for m = 1:M
    M_tf = [R(m,3).*out(m,1:3) R(m,2).*out(m,4:6)] + M_tf;
end

%% Dispaly Results
disp('Synchrony results');
disp(out); disp(M_tf);
