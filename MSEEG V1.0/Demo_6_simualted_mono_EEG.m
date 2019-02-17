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

%% Parameters
N  = 768;                 % Total number of time samples
M  = 1024;                % Total number of frequency bins
fs = 32;                  % Sampling frequency
t  = 0:1/fs:(N-1)/fs;     % Time array
f  = 0:fs/(2*M-1):fs/2;   % Positive frequency array
rng(34);                  % ensuring reproducibility

%% Creating mono-channel waveforms
Xb = EEG_back(N, fs);
Xs = EEG_seiz(N, fs);

%% Frequency Analysis
Xbf = abs(fftshift(fft(Xb, 2*M))); Xbf = Xbf(M+1:2*M);
Xsf = abs(fftshift(fft(Xs, 2*M))); Xsf = Xsf(M+1:2*M);

%% Time-Frequency Analysis
Xbtf = abs(quadtfd(Xb, N/4-1, 1, 'emb',0.075,0.5,M));
Xstf = abs(quadtfd(Xs, N/4-1, 1, 'emb',0.025,0.5,M));

%% Plotting
K = 7; cnt = 0; tfr = 6;
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

figure('Color',[1,1,1],'Position',[100 0 800 650]);  colormap(1-gray);
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
    print(1,'sim_mono_back','-dpdf','-r512');
    print(2,'sim_mono_seiz','-dpdf','-r512');
    print(3,'sim_mono_back_flat','-dpdf','-r512');
    print(4,'sim_mono_seiz_flat','-dpdf','-r512');
else
    return
end
