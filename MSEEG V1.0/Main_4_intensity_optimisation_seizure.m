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
addpath(genpath('Data\Power Maps'));
rng(1);

%% Parameters
M  = 7;     % Total number of patients
n  = 1;     % Number of source signals
K  = 100;   % Number of initial intensities

%% Clerc & Kennedy Constriction Coefficients
kappa = 1;             % 0 <= kappa <=1 (defualt 1)
phi1  = 2.05;
phi2  = 2.05;
phi   = phi1 + phi2;   % phi >= 4 (default 4.1)
chi   = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

%% Parameters of PSO
params.MaxIt = 50;           % Maximum number of iterations
params.nPop  = 30;           % Population size (Swarm size)
params.w     = chi;          % Inertia coefficient
params.wdamp = 0.95;         % Damping ratio of inertia coefficient
params.c1    = chi*phi1;     % Personal acceleration coefficient
params.c2    = chi*phi2;     % Social (Global) acceleration coefficient
params.ShowIterInfo = 1;     % Flag for displaying iteration information

%% Calculating Power Threshold
disp('Calculating Power Threshold')
P = [];
for m = 1:M
    % Loading Real EEG
    load(['Patient_' num2str(m) '_Data_corrected.mat'],'EEG_S');
    % Multi-Channel Power
    [Ns, ch_n, N] = size(EEG_S);
    Ps = zeros(ch_n, Ns);
    for i = 1:Ns
        tmp(:,:) = EEG_S(i,:,:);
        Ps(:,i)  = sum(tmp.^2,2)./N;
    end
    P = [P max(Ps)];
end

%% Defining maximum power of interest
mP = linspace(min(P),max(P),K);

%% Optimisation
for m = 1:M
    disp(['Patient No. ' num2str(m)])
    % Loading
    load(['Patient_' num2str(m) '_power.mat'],'P_s');
    P = P_s;
    % Initialise
    [ch_n, Ns] = size(P);
    Error = zeros(1,K);
    Corr  = zeros(1,K);
    % Optimise
    for k = 1:K
        Lamda_s = zeros(ch_n, Ns);
        for ns = 1:Ns
            out = PSO(params, n, P(:,ns), mP(k));
            opt_pos_seiz  = out.BestSol.Position;
            Lamda_s(:,ns) = multi_prop_opt(opt_pos_seiz, mP(k));
        end
        Error(1,k) = 100*norm(P - Lamda_s)/norm(P);
        Corr(1,k)  = 100*abs(corr(P(:),Lamda_s(:)));
        disp(100*k/K)
    end
    % Saving
    save(['Data/Optimisation/Intensity/Patient_' num2str(m) '_seizure'],'Error','Corr','mP');
end