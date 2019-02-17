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
% NOTE: This script requires the database of real multi-sensor newborn EEG.

%% Initialisation
clear; close all; clc; warning off;
addpath(genpath('Toolbox'));
addpath(genpath('Data\EEG'));

%% Loading EEG
load('Patient_6_Data_corrected.mat','EEG_B','EEG_S','fs','Channel_label');

%% Selecting one multi-channel Segment
B_seg = 94:96;
S_seg = 2:4;
back = []; seiz = [];
for i = 1:length(S_seg)
    tmp_B(:,:) = EEG_B(B_seg(i),:,:);
    tmp_S(:,:) = EEG_S(S_seg(i),:,:);
    back = [back tmp_B];
    seiz = [seiz tmp_S];
end

%% Parameters
[ch_n, N] = size(back);
M  = 1024;
t  = 0:1/fs:(N-1)/fs;
f  = 0:fs/(2*M-1):fs/2;

%% Creating mono-channel waveforms
Xb = mean(back);
Xs = mean(seiz);

%% Frequency Analysis
Xbf = abs(fftshift(fft(Xb, 2*M))); Xbf = Xbf(M+1:2*M);
Xsf = abs(fftshift(fft(Xs, 2*M))); Xsf = Xsf(M+1:2*M);

%% Time-Frequency Analysis
Xbtf = abs(quadtfd(Xb, N/4-1, 1, 'emb',0.075,0.5,M));
Xstf = abs(quadtfd(Xs, N/4-1, 1, 'emb',0.025,0.5,M));

%% Plotting
K = 7; cnt = 0; tfr = 6; sh = 40;
figure('Color',[1,1,1],'Position',[100 -50 900 700]);
plot_multichannel(back, sh, fs,[0 0.447058826684952 0.74117648601532],2); hold on

plot(t, back(7,:) - sh.*(7-1),'color','k','linewidth',2); hold on
plot(t, back(10,:) - sh.*(10-1),'color','k','linewidth',2); hold on
plot(t, back(8,:) - sh.*(8-1),'color','r','linewidth',2); hold on
plot(t, back(11,:) - sh.*(11-1),'color','r','linewidth',2); hold on
plot(t, back(13,:) - sh.*(13-1),'color',[0.6 0.6 0],'linewidth',2); hold on
plot(t, back(16,:) - sh.*(16-1),'color',[0.6 0.6 0],'linewidth',2); hold on
plot(t, back(14,:) - sh.*(14-1),'color','m','linewidth',2); hold on
plot(t, back(15,:) - sh.*(15-1),'color','m','linewidth',2); hold on

y = [(sh*(1-ch_n))-sh sh sh (sh*(1-ch_n))-sh];
x = [t(1) t(1) t(1.5*fs) t(1.5*fs)];
h = fill(x, y, 'g'); set(h,'FaceAlpha',0.2,'Linestyle','none'); hold on;
x = [t(3.5*fs) t(3.5*fs) t(4.5*fs) t(4.5*fs)];
h = fill(x, y, 'g'); set(h,'FaceAlpha',0.2,'Linestyle','none');
x = [t(8*fs) t(8*fs) t(24*fs) t(24*fs)];
h = fill(x, y, 'g'); set(h,'FaceAlpha',0.2,'Linestyle','none');
hold off; xlabel('Time (s)'); xlim([0 t(end)]);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
set(gca,'yticklabel',flipud(Channel_label),'fontweight','bold','fontsize',18,'FontName','Times');

figure('Color',[1,1,1],'Position',[100 -50 900 700]);
plot_multichannel(seiz, sh, fs,'r',2); xlabel('Time (s)'); hold on;
plot(t, seiz(3,:) - sh.*(2),'color','k','linewidth',2); hold on
plot(t, seiz(9,:) - sh.*(8),'color','k','linewidth',2); hold on
plot(t, seiz(14,:) - sh.*(13),'color','k','linewidth',2); hold on
plot(t, seiz(18,:) - sh.*(17),'color','k','linewidth',2); hold off
set(gcf,'Units','inches'); screenposition = get(gcf,'Position'); xlim([0 t(end)]);
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
set(gca,'yticklabel',flipud(Channel_label),'fontweight','bold','fontsize',18,'FontName','Times');

kk = zeros(K,K);
for a = 1:K
    for b = 1:K
        cnt = cnt + 1;
        kk(a,b) = cnt;
    end
end
figure('Color',[1,1,1],'Position',[100 0 800 650]); colormap(1-gray);
subplot(K,K,kk(1:end-1,1)); plot(Xb,t,'linewidth',1.5);
grid on; axis([min(Xb)-0.5 max(Xb)+0.5 0 t(end)]);
set(gca,'Xdir', 'reverse','FontWeight','bold','FontSize',18,'xticklabel',{})
ylabel('Time (s)','FontSize',18);
subplot(K,K,kk(end,2:end));
plot(f,Xbf,'r','linewidth',1.5);
grid on; axis([0 fs/2 0 max(Xbf)+10]);
set(gca,'FontWeight','bold','FontSize',18,'yticklabel',{})
xlabel('Frequency (Hz)','FontSize',18);
tmp = kk(1:end-1,2:end);
subplot(K,K,tmp(:));
imagesc(f,t(1:tfr:end),abs(Xbtf(:,1:tfr:end))'); axis xy;
set(gca,'FontWeight','bold','xticklabel',{},'yticklabel',{})
title('Joint Time-Frequency Domain','FontSize',20)
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[100 0 800 650]); colormap(1-gray);
subplot(K,K,kk(1:end-1,1)); plot(Xs,t,'linewidth',1.5);
grid on; axis([min(Xs)-0.5 max(Xs)+0.5 0 t(end)]);
set(gca,'Xdir', 'reverse','FontWeight','bold','FontSize',18,'xticklabel',{})
ylabel('Time (s)','FontSize',18);
subplot(K,K,kk(end,2:end));
plot(f,Xsf,'r','linewidth',1.5);
grid on; axis([0 fs/2 0 max(Xsf)+10]);
set(gca,'FontWeight','bold','FontSize',18,'yticklabel',{})
xlabel('Frequency (Hz)','FontSize',18);
tmp = kk(1:end-1,2:end);
subplot(K,K,tmp(:));
imagesc(f,t(1:tfr:end),abs(Xstf(:,1:tfr:end))'); axis xy;
set(gca,'FontWeight','bold','xticklabel',{},'yticklabel',{})
title('Joint Time-Frequency Domain','FontSize',20)
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[100 0 800 650]);
subplot(K,K,kk(1:end-1,1)); plot(Xb,t,'linewidth',1.5);
grid on; axis([min(Xb)-0.5 max(Xb)+0.5 0 t(end)]);
set(gca,'Xdir', 'reverse','FontWeight','bold','FontSize',18,'xticklabel',{})
ylabel('Time (s)','FontSize',18);
subplot(K,K,kk(end,2:end));
plot(f,Xbf,'r','linewidth',1.5);
grid on; axis([0 fs/2 0 max(Xbf)+10]);
set(gca,'FontWeight','bold','FontSize',18,'yticklabel',{})
xlabel('Frequency (Hz)','FontSize',18);
tmp = kk(1:end-1,2:end);
subplot(K,K,tmp(:));
flatwf(f,t(1:tfr:end),abs(Xbtf(:,1:tfr:end))','w','k');
set(gca,'FontWeight','bold','xticklabel',{},'yticklabel',{})
title('Joint Time-Frequency Domain','FontSize',20)
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[100 0 800 650]);
subplot(K,K,kk(1:end-1,1)); plot(Xs,t,'linewidth',1.5);
grid on; axis([min(Xs)-0.5 max(Xs)+0.5 0 t(end)]);
set(gca,'Xdir', 'reverse','FontWeight','bold','FontSize',18,'xticklabel',{})
ylabel('Time (s)','FontSize',18);
subplot(K,K,kk(end,2:end));
plot(f,Xsf,'r','linewidth',1.5);
grid on; axis([0 fs/2 0 max(Xsf)+10]);
set(gca,'FontWeight','bold','FontSize',18,'yticklabel',{})
xlabel('Frequency (Hz)','FontSize',18);
tmp = kk(1:end-1,2:end);
subplot(K,K,tmp(:));
flatwf(f,t(1:tfr:end),abs(Xstf(:,1:tfr:end))','w','k');
set(gca,'FontWeight','bold','xticklabel',{},'yticklabel',{})
title('Joint Time-Frequency Domain','FontSize',20)
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'real_multi_back','-dpdf','-r512');
    print(2,'real_multi_seiz','-dpdf','-r512');
    print(3,'real_mono_back','-dpdf','-r512');
    print(4,'real_mono_seiz','-dpdf','-r512');
    print(5,'real_mono_back_flat','-dpdf','-r512');
    print(6,'real_mono_seiz_flat','-dpdf','-r512');
else
    return
end
