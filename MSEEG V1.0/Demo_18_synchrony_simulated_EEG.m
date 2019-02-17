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
clear; close all; clc; warning off
addpath(genpath('Toolbox'));
addpath(genpath('Data\Synchrony'));
addpath(genpath('Data\EEG'));

%% Parameters
per = 0.80;   % Causality plotting threshold
M   = 7;      % Total number of patients

%% Averaging correlation matrices
load('Ratios');
Synchrony_B = 0;
Synchrony_S = 0;
for m = 1:M
    disp(['Patient No. ' num2str(m)]);
    % Loading Real EEG
    load(['simulated_EEG_' num2str(m) '_synch.mat'],'CB','CS','Channel_label');
    % Averaging
    disp('Averaging all patients ...');
    Synch_B = 0;
    Nb = length(CB);
    for nb = 1:Nb
        Synch_B  = Synch_B  + abs(CB{1,nb});
    end
    Synchrony_B  = R(m,3)*(Synch_B./Nb  + Synchrony_B);
    
    Synch_S = 0;
    Ns = length(CS);
    for ns = 1:Ns
        Synch_S  = Synch_S  + abs(CS{1,ns});
    end
    Synchrony_S  = R(m,2)*(Synch_S./Ns  + Synchrony_S);
end
Synchrony_B  = Synchrony_B./max(diag(Synchrony_B));
Synchrony_S  = Synchrony_S./max(diag(Synchrony_S));

%% Thresholding
disp('Thresholding ...');
Synch_B_thr  = Synchrony_B.*(Synchrony_B>(1-per));
Synch_S_thr  = Synchrony_S.*(Synchrony_S>(1-per));

%% Plotting
figure('Color',[1 1 1],'Position',[100, 50, 550 450]);
colormap(1-gray); imagesc(abs(Synch_B_thr)); axis xy; grid on; caxis([0 1]);
set(gca,'XTickLabel',Channel_label,'YTickLabel',Channel_label);
set(gca,'Xtick',1:20,'FontWeight','bold','FontSize',12,'Ytick',1:20,'FontWeight','bold','FontSize',12);
title('Background averaged synchrony matrix','FontSize',16);
rotateXLabels(gca(),90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100, 50, 550 450]);
colormap(1-gray); imagesc(abs(Synch_S_thr)); axis xy; grid on; caxis([0 1]);
set(gca,'XTickLabel',Channel_label,'YTickLabel',Channel_label);
set(gca,'Xtick',1:20,'FontWeight','bold','FontSize',12,'Ytick',1:20,'FontWeight','bold','FontSize',12);
title('Seizure averaged synchrony matrix','FontSize',16);
rotateXLabels(gca(),90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'sim_synchrony_back','-dpdf','-r400');
    print(2,'sim_synchrony_seiz','-dpdf','-r400');
else
    return
end
