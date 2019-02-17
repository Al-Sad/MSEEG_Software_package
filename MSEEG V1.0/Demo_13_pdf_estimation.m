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
addpath(genpath('Data\Optimisation\Noise'));

%% Loading Results
load('Noise_Estimation');
load('Clean_PDF');
load('Noisy_PDF');

%% PD Estimation
x = -70:0.1:70;
P_eeg_s = pdf('tLocationScale',x,m_s,s_s,n_s);
P_sim_s = pdf('tLocationScale',x,m_sim_s,s_sim_s,n_sim_s);
P_nos_s = pdf('tLocationScale',x,m_nos_s,s_nos_s,n_nos_s);
P_eeg_b = pdf('tLocationScale',x,m_b,s_b,n_b);
P_sim_b = pdf('tLocationScale',x,m_sim_b,s_sim_b,n_sim_b);
P_nos_b = pdf('tLocationScale',x,m_nos_b,s_nos_b,n_nos_b);

%% Display
disp('Real seizures');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_s,s_s,n_s);
disp('Simualted clean seizures');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_sim_s,s_sim_s,n_sim_s);
disp('Simualted noisy seizures');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_nos_s,s_nos_s,n_nos_s);

disp('Real background');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_b,s_b,n_b);
disp('Simualted clean background');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_sim_b,s_sim_b,n_sim_b);
disp('Simualted noisy background');
fprintf('m = %0.2f, s = %0.2f, v = %0.2f\n',m_nos_b,s_nos_b,n_nos_b);

%% Plotting
figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_s,H_s,1); axis([-70 70 0 0.05]); hold on;
plot(x,P_eeg_s,'r','linewidth',3); hold on;
plot(repmat(m_s,1,100), linspace(0,0.05,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_sim_s,H_sim_s,1); axis([-70 70 0 0.05]); hold on;
plot(x,P_sim_s,'r','linewidth',3); hold on;
plot(repmat(m_sim_s,1,100), linspace(0,0.05,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_nos_s,H_nos_s,1); axis([-70 70 0 0.05]); hold on;
plot(x,P_nos_s,'r','linewidth',3); hold on;
plot(repmat(m_nos_s,1,100), linspace(0,0.05,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_b,H_b,1); axis([-70 70 0 0.12]); hold on;
plot(x,P_eeg_b,'r','linewidth',3); hold on;
plot(repmat(m_b,1,100), linspace(0,0.12,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_sim_b,H_sim_b,1); axis([-70 70 0 0.12]); hold on;
plot(x,P_sim_b,'r','linewidth',3); hold on;
plot(repmat(m_sim_b,1,100), linspace(0,0.12,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(x_nos_b,H_nos_b,1); axis([-70 70 0 0.12]); hold on;
plot(x,P_nos_b,'r','linewidth',3); hold on;
plot(repmat(m_nos_b,1,100), linspace(0,0.12,100),'r--','linewidth',2); grid on;
legend('Data histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Amplitude','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'real_seizure_PDF','-dpdf','-r400');
    print(2,'simulated_clean_seizure_PDF','-dpdf','-r400');
    print(3,'simulated_noisy_seizure_PDF','-dpdf','-r400');
    print(4,'real_background_PDF','-dpdf','-r400');
    print(5,'simulated_clean_background_PDF','-dpdf','-r400');
    print(6,'simulated_noisy_background_PDF','-dpdf','-r400');
end
