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
clear; close all; clc; warning off;
addpath(genpath('Toolbox'));
addpath(genpath('Data\Correlation'));
addpath(genpath('Data\EEG'));

%% Parameters
per  = 0.5;   % Correlation plotting threshold
M    = 7;     % Total number of patients

%% Averaging correlation matrices
load('Ratios');
Correlation_B_t = 0; Correlation_B_tf = 0;
Correlation_S_t = 0; Correlation_S_tf = 0;
for m = 1:M
    disp(['Patient No. ' num2str(m)]);
    % Loading Real EEG
    load(['Patient_' num2str(m) '_corr.mat'],'CB_t','CS_t','CB_tf','CS_tf','Channel_label');
    % Averaging
    disp('Averaging all patients ...');
    Corr_B_t = 0; Corr_B_tf = 0; 
    Nb = length(CB_t);
    for nb = 1:Nb
        Corr_B_t  = Corr_B_t  + CB_t{1,nb};
        Corr_B_tf = Corr_B_tf + CB_tf{1,nb};
    end
    Correlation_B_t  = R(m,3)*(Corr_B_t./Nb  + Correlation_B_t);
    Correlation_B_tf = R(m,3)*(Corr_B_tf./Nb + Correlation_B_tf);
    
    Corr_S_t = 0; Corr_S_tf = 0;
    Ns = length(CS_t);
    for ns = 1:Ns
        Corr_S_t  = Corr_S_t  + CS_t{1,ns};
        Corr_S_tf = Corr_S_tf + CS_tf{1,ns};
    end
    Correlation_S_t  = R(m,2)*(Corr_S_t./Ns  + Correlation_S_t);
    Correlation_S_tf = R(m,2)*(Corr_S_tf./Ns + Correlation_S_tf); 
end
Correlation_B_t  = Correlation_B_t./max(diag(Correlation_B_t));
Correlation_B_tf = Correlation_B_tf./max(diag(Correlation_B_tf));
Correlation_S_t  = Correlation_S_t./max(diag(Correlation_S_t));
Correlation_S_tf = Correlation_S_tf./max(diag(Correlation_S_tf));

%% Thresholding
disp('Thresholding ...');
Corr_B_t_thr  = Correlation_B_t.*(Correlation_B_t>(1-per));
Corr_B_tf_thr = Correlation_B_tf.*(Correlation_B_tf>(1-per));
Corr_S_t_thr  = Correlation_S_t.*(Correlation_S_t>(1-per));
Corr_S_tf_thr = Correlation_S_tf.*(Correlation_S_tf>(1-per));

%% Plotting
figure('Color',[1 1 1],'Position',[100, 50, 550 450]);
colormap(1-gray); imagesc(Corr_B_tf_thr); axis xy; grid on; caxis([0 1]);
set(gca,'XTickLabel',Channel_label,'YTickLabel',Channel_label);
set(gca,'Xtick',1:20,'FontWeight','bold','FontSize',12,'Ytick',1:20,'FontWeight','bold','FontSize',12);
title('Background averaged TF correlation matrix','FontSize',16);
rotateXLabels(gca(),90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1 1 1],'Position',[100, 50, 550 450]);
colormap(1-gray); imagesc(Corr_S_tf_thr); axis xy; grid on; caxis([0 1]);
set(gca,'XTickLabel',Channel_label,'YTickLabel',Channel_label);
set(gca,'Xtick',1:20,'FontWeight','bold','FontSize',12,'Ytick',1:20,'FontWeight','bold','FontSize',12);
title('Seizure averaged TF correlation matrix','FontSize',16);
rotateXLabels(gca(),90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'real_corr_tf_back','-dpdf','-r400');
    print(2,'real_corr_tf_seiz','-dpdf','-r400');
else
    return
end
