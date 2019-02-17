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
addpath(genpath('Data\Optimisation\Finalised'));

%% Parameters
M  = 7;    % Total number of patients
sh = 200;  % Multi-sensor EEG plotting vertical spacing

%% Main
for m = 1:M
    % Loading Data
    disp(['Patient # ' num2str(m)]);
    load(['Data\Simulated EEG\Patient_' num2str(m)]);
    t = 0:1/fs:(length(X_eeg)-1)/fs;
    % Plotting
    figure('Color',[1,1,1],'Position',[0 -200 1000 850]);
    subplot(20,1,1:18)
    plot_multichannel(X_sim_noisy, sh, fs,'r'); hold on; xlim([t(1) t(end)]);
    set(gca,'yticklabel',Channel_label,'xticklabel','','fontweight','bold','fontsize',20,'FontName','Times');
%     title('Simulated multi-sensor neonatal EEG','fontsize',30);
    subplot(20,1,19:20)
    plot(t,mask,'k','linewidth',2); grid on; axis([t(1) t(end) -0.25 1.25]);
    set(gca,'Ytick',0:1,'yticklabel',{'B','S'},'FontWeight','bold','FontSize',20,'FontName','Times');
    xlabel('Time (s)','fontweight','bold','fontsize',20);
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
    % Saving
    opt = input('Do you want to save results (Y/N)\n','s');
    if(opt == 'y' || opt == 'Y')
        print(1,['simulated_multi_EEG_' num2str(m)],'-dpdf','-r400');
        close all;
    else
        close all;
    end
end