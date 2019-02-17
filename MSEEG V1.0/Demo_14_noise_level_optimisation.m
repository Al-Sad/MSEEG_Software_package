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

%% Parameters
bin = 16;   % Total number of bins in a histogram

%% Loading Results
load('Noise_Level');

%% Estimating noise level PDF
L_s = fitdist(N_opt_s','ExtremeValue');
L_b = fitdist(N_opt_b','InverseGaussian');
x = 0:0.1:100;
Gs = pdf('ExtremeValue',x,L_s.mu,L_s.sigma);
Gb = pdf('InverseGaussian',x,L_b.mu,L_b.lambda);
[Hs, Xs] = hist(N_opt_s,bin);
Hs = Hs ./ trapz(Xs,Hs);
[Hb, Xb] = hist(N_opt_b,bin);
Hb = Hb ./ trapz(Xb,Hb);

%% Saving
D_s = 'ExtremeValue';
D_b = 'InverseGaussian';
M_s = L_s.mu;
S_s = L_s.sigma;
M_b = L_b.mu;
S_b = L_b.lambda;
save('Data\Optimisation\Noise\Noise_Level_distribution','D_s','D_b',...
    'M_s','S_s','M_b','S_b');

%% Display
disp('Seizure noise level')
fprintf('m = %0.2f, s = %0.2f\n',M_s,S_s);
disp('Background noise level')
fprintf('m = %0.2f, s = %0.2f\n',M_b,S_b);

%% Plotting
figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(Xs,Hs,1); axis([0 100 0 0.03]); hold on;
plot(x,Gs,'r','linewidth',3); hold on;
plot(repmat(L_s.mu,1,100), linspace(0,0.04,100),'r--','linewidth',2); grid on;
legend('Noise level histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Noise Level (%)','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1],'Position',[0 0 650 550]);
bar(Xb,Hb,1); axis([0 100 0 0.35]); hold on;
plot(x,Gb,'r','linewidth',3); hold on;
plot(repmat(L_b.mu,1,100), linspace(0,0.36,100),'r--','linewidth',2); grid on;
legend('Noise level histogram','Estimated PDF')
set(gca,'fontweight','bold','fontsize',16,'FontName','Times');
xlabel('Noise Level (%)','fontsize',18); ylabel('Probability Density','fontsize',16);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save all results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'noise_level_seizure','-dpdf','-r400');
    print(2,'noise_level_background','-dpdf','-r400');
end