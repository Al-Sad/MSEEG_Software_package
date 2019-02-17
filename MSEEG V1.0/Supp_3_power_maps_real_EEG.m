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

%% Parameters
M = 7;   % Total number of patients

%% Main
for m = 1:M
    % Loading original power maps
    load(['Power Maps\Patient_' num2str(m) '_power.mat'],'P_eeg','mask_seg','Channel_label');
    [ch_n, N] = size(P_eeg);
    % Plotting
    figure('Color',[1,1,1],'Position',[100 0 700 700]);
    subplot(10,1,1:9)
    imagesc(10*log10(P_eeg)); axis xy; colormap jet; xlim([1 N]);
    g = colorbar('location','NorthOutside'); ylabel(g, 'Power (dB)','fontweight','bold','fontsize',18);
    caxis([min(min(10*log10(P_eeg))) max(max(10*log10(P_eeg)))])
    set(gca,'yticklabel',Channel_label,'Ytick',1:20,'xticklabel','','fontweight','bold','fontsize',18,'FontName','Times');
%     title(['Patient ' num2str(m) ' power map'],'FontWeight','bold','FontSize',18);
    subplot(10,1,10)
    plot(mask_seg,'k','linewidth',2); grid on;
    set(gca,'Ytick',0:1,'yticklabel',{'B','S'},'FontWeight','bold','FontSize',16);
    xlabel('Segment Number','fontweight','bold','fontsize',18); axis([1 N -0.25 1.25]);
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
end

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    for m = 1:M
        print(m,['real_power_map_' num2str(m)],'-dpdf','-r400');
    end
end