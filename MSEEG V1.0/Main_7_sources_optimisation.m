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
addpath(genpath('Data\Optimisation\Sources'));
addpath(genpath('Data\Power Maps'));
addpath(genpath('Data\EEG'));

%% Parameters
M   = 7;                      % Total number of patients
N   = 6;                      % Maximum number of source signals
Pop = [30 50 75 100 150 200]; % Population size (Swarm size)

%% Loading optimal intensities
load('Data\Optimisation\Intensity\Optimal_Intensities');

%% Computing Correlations and Errors
disp('Computing Correlations and Errors')
Error_s = cell(1,M);
Error_b = cell(1,M);
Corr_s  = cell(1,M);
Corr_b  = cell(1,M);
for m = 1:M
    % Loading real power maps
    disp(['Patient # ' num2str(m)]);
    load(['Patient_' num2str(m) '_power.mat'],'P_s','P_b');
    [ch_n, ~] = size(P_s);
    for p = 1:length(Pop)
        if(p == 5 && m == 4)
            N = 5;
        elseif(p == 6 && m == 4)
            N = 5;
        else
            N = 6;
        end
        for n = 1:N
            % Loading Results
            load(['Patient_' num2str(m) '_Pop_' num2str(Pop(p)) '_n' num2str(n) '.mat'],...
                'opt_pos_back','opt_pos_seiz');
            % Computing simulated power maps
            Ns = length(opt_pos_seiz);
            Nb = length(opt_pos_back);
            Lamda_s = zeros(ch_n, Ns);
            Lamda_b = zeros(ch_n, Nb);
            for ns = 1:Ns
                Xs = opt_pos_seiz{1,ns}(:,1);
                Ys = opt_pos_seiz{1,ns}(:,2);
                Zs = opt_pos_seiz{1,ns}(:,3);
                Lamda_s(:,ns) = multi_prop_opt([Xs Ys Zs],mP_s(m));
            end
            for nb = 1:Nb
                Xb = opt_pos_back{1,nb}(:,1);
                Yb = opt_pos_back{1,nb}(:,2);
                Zb = opt_pos_back{1,nb}(:,3);
                Lamda_b(:,nb) = multi_prop_opt([Xb Yb Zb],mP_b(m));
            end
            Error_s{1,m}(n,p) = 100*norm(P_s - Lamda_s)/norm(P_s);
            Corr_s{1,m}(n,p)  = 100*abs(corr(reshape(P_s,1,Ns*ch_n)',reshape(Lamda_s,1,Ns*ch_n)'));
            Error_b{1,m}(n,p) = 100*norm(P_b - Lamda_b)/norm(P_b);
            Corr_b{1,m}(n,p)  = 100*abs(corr(reshape(P_b,1,Nb*ch_n)',reshape(Lamda_b,1,Nb*ch_n)'));
        end
    end
end

%% Computing and saving optimal results
p  = zeros(M,1);
Np = zeros(M,1);
E  = zeros(M,1);
C  = zeros(M,1);
for m = 1:M
    PP = 0.5.*(100 + Corr_b{1,m} - Error_b{1,m});
    [Vr, Ir] = max(PP);
    [Vc, Ic] = max(Vr);
    p(m,1)   = Ir(Ic);
    Np(m,1)  = Pop(Ic);
    E(m,1)   = Error_b{1,m}(Ir(Ic), Ic);
    C(m,1)   = Corr_b{1,m}(Ir(Ic), Ic);
    load(['Patient_' num2str(m) '_Pop_' num2str(Np(m,1)) '_n' num2str(p(m,1)) '.mat'],'opt_pos_back');
    save(['Data\Optimisation\Positions\Optimal Positions\Patient_' num2str(m)],'opt_pos_back');
end
save('Data\Optimisation\Positions\Source_Optimisation_Background','p','Np','E','C');

p  = zeros(M,1);
Np = zeros(M,1);
E  = zeros(M,1);
C  = zeros(M,1);
for m = 1:M
    PP = 0.5.*(100 + Corr_s{1,m} - Error_s{1,m});
    [Vr, Ir] = max(PP);
    [Vc, Ic] = max(Vr);
    p(m,1)   = Ir(Ic);
    Np(m,1)  = Pop(Ic);
    E(m,1)   = Error_s{1,m}(Ir(Ic), Ic);
    C(m,1)   = Corr_s{1,m}(Ir(Ic), Ic);
    load(['Patient_' num2str(m) '_Pop_' num2str(Np(m,1)) '_n' num2str(p(m,1)) '.mat'],'opt_pos_seiz');
    save(['Data\Optimisation\Positions\Optimal Positions\Patient_' num2str(m)],'opt_pos_seiz','-append');
end
save('Data\Optimisation\Positions\Source_Optimisation_Seizure','p','Np','E','C');
