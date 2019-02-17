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
addpath(genpath('Data\Power Maps'));
rng(1);

%% Parameters
M    = 7;                      % Total number of patients
N    = 6;                      % Maximum number of source signals
Pop  = [30 50 75 100 150 200]; % Population size (Swarm size)

%% Clerc & Kennedy Constriction Coefficients
kappa = 1;             % 0 <= kappa <=1 (defualt 1)
phi1  = 2.05;
phi2  = 2.05;
phi   = phi1 + phi2;   % phi >= 4 (default 4.1)
chi   = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));

%% Parameters of PSO
params.MaxIt = 300;            % Maximum number of iterations
params.w     = chi;            % Inertia coefficient
params.wdamp = 0.95;           % Damping ratio of inertia coefficient
params.c1    = chi*phi1;       % Personal acceleration coefficient
params.c2    = chi*phi2;       % Social (Global) acceleration coefficient

%% Loading optimal intensities
load('Data\Optimisation\Intensity\Optimal_Intensities');

%% Optimisation
for m = 1:M
    disp(['Patient No. ' num2str(m)])
    load(['Patient_' num2str(m) '_power.mat']);
    [ch_n, Ns] = size(P_s);
    [~, Nb]    = size(P_b);
    for i = 1:length(Pop)
        params.nPop = Pop(i);
        disp(['Number of populations ' num2str(Pop(i))])
        for n = 1:N
            opt_pos_seiz  = cell(1,Ns);
            opt_cost_seiz = zeros(1,Ns);
            opt_pos_back  = cell(1,Nb);
            opt_cost_back = zeros(1,Nb);
            disp(['Number of sources ' num2str(n)])
            disp('Optimising seizure segments');
            for ns = 1:Ns
                out = PSO(params, n, P_s(:,ns), mP_s(m));
                opt_pos_seiz{1,ns}  = out.BestSol.Position;
                opt_cost_seiz(1,ns) = out.BestSol.Cost;
                disp(100*ns/Ns)
            end
            disp('Optimising background segments');
            for nb = 1:Nb
                out = PSO(params, n, P_b(:,nb), mP_b(m));
                opt_pos_back{1,nb}  = out.BestSol.Position;
                opt_cost_back(1,nb) = out.BestSol.Cost;
                disp(100*nb/Nb)
            end
            % Saving
            disp('Saving ...')
            filename = ['Data\Optimisation\Sources\Population ' num2str(Pop(i)) '\Patient_' num2str(m) ...
                '_Pop_' num2str(Pop(i)) '_n' num2str(n)];
            save(filename,'opt_pos_seiz','opt_pos_back','opt_cost_seiz','opt_cost_back');
        end
    end
end
