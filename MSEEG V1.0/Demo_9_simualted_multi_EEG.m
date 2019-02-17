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

%% Parameters
N     = 2048;             % Total number of samples
fs    = 32;               % Sampling frequency
t     = 0:1/fs:(N-1)/fs;  % Time array
R     = 5.95;             % Head Radius in cm
mP_b  = 100;              % Initial intensity of background EEG
mP_s  = 300;              % Initial intensity of seizure EEG
speed = 2;                % Propagation speed (cm/s)
sh    = 25;               % Multi-sensor plotting vertical spacing
rng(3);                   % Ensuring reproducibility

%% Loading Channel Labels
load('Channels');

%% EEG Mask
mask = [zeros(1,N/8) ones(1,5*N/8) zeros(1,2*N/8)];

%% Originating Locations
event_b{1} = [0.2 -0.4 -0.6 ; 0.8 0.4 -1.0 ; -0.6 0.4 0.5];
event_s{1} = [-3 -2.5 1.5 ; -2 1 3];

%% Multi-sensor EEG Generation
[X, Xb, Xs] = multi_sensor_EEG(mask, fs, event_b, event_s, mP_b, mP_s, speed);
[ch_n, ~] = size(X);

%% Plotting Time-Space Analysis
figure('Color',[1,1,1],'Position',[0 -200 1000 850]);
plot_multichannel(Xb, sh, fs,'b'); hold on; xlim([t(1) t(end)]);
set(gca,'yticklabel',Channel_label,'fontweight','bold','fontsize',22,'FontName','Times');
xlabel('Time (s)','fontweight','bold','fontsize',22);
title('Multi-sensor synthetic neonatal background','fontsize',30);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 -200 1000 850]);
plot_multichannel(Xs, sh, fs,'r'); hold on; xlim([t(1) t(end)]);
set(gca,'yticklabel',Channel_label,'fontweight','bold','fontsize',22,'FontName','Times');
xlabel('Time (s)','fontweight','bold','fontsize',22);
title('Multi-sensor synthetic neonatal seizure','fontsize',30);
y = [(sh*(1-ch_n))-sh sh sh (sh*(1-ch_n))-sh];
if(sum(mask)==N), x = [0 0 t(end) t(end)]; h = fill(x, y, 'k'); set(h,'FaceAlpha',0.1,'Linestyle','none');
else
    cnt1 = 1; cnt2 = 1;
    temp = [];
    flag1 = 0;
    for i = 1:N
        if(mask(i) == 1)
            temp(1,cnt2) = i;
            cnt2 = cnt2 + 1;
            flag1 = 1;
        elseif(mask(i) == 0 && flag1)
            ind{cnt1} = temp;
            temp = [];
            cnt1 = cnt1 + 1;
            cnt2 = 1;
            flag1 = 0;
        end
    end
    for i = 1:cnt1-1
        x = [t(ind{i}(1)) t(ind{i}(1)) t(ind{i}(end)) t(ind{i}(end))];
        h = fill(x, y, 'g'); set(h,'FaceAlpha',0.2,'Linestyle','none')
    end
end
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 -200 1000 850]);
plot_multichannel(X, sh, fs,'b'); hold on; xlim([t(1) t(end)]);
set(gca,'yticklabel',Channel_label,'fontweight','bold','fontsize',22,'FontName','Times');
xlabel('Time (s)','fontweight','bold','fontsize',22);
title('Multi-sensor synthetic neonatal EEG','fontsize',30)
y = [(sh*(1-ch_n))-sh sh sh (sh*(1-ch_n))-sh];
if(sum(mask)==N), x = [0 0 t(end) t(end)]; h = fill(x, y, 'k'); set(h,'FaceAlpha',0.1,'Linestyle','none');
else
    cnt1 = 1; cnt2 = 1;
    temp = [];
    flag1 = 0;
    for i = 1:N
        if(mask(i) == 1)
            temp(1,cnt2) = i;
            cnt2 = cnt2 + 1;
            flag1 = 1;
        elseif(mask(i) == 0 && flag1)
            ind{cnt1} = temp;
            temp = [];
            cnt1 = cnt1 + 1;
            cnt2 = 1;
            flag1 = 0;
        end
    end
    for i = 1:cnt1-1
        x = [t(ind{i}(1)) t(ind{i}(1)) t(ind{i}(end)) t(ind{i}(end))];
        h = fill(x, y, 'g'); set(h,'FaceAlpha',0.2,'Linestyle','none')
    end
end
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'sim_multi_back','-dpdf','-r512');
    print(2,'sim_multi_seiz','-dpdf','-r512');
    print(3,'sim_multi_EEG','-dpdf','-r512');
else
    return
end