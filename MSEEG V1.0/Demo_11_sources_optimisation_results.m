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
addpath(genpath('Data\Optimisation\Positions'));
addpath(genpath('Data\EEG'));

%% Parameters
M = 7;             % Total number of patients
brain_R = 4.76;    % Brain Radius in cm

%% Loading Database information
load('Ratios');

%% Display Results
load('Source_Optimisation_Background');
fprintf('\t\t\t\t\tBackground Results\n');
fprintf('\t  #\t\t\tp\t\tN_p\t\t  Corr\t\t Error\n');
for m = 1:M
    fprintf('\t  %1.0f\t\t\t%1.0f\t\t%1.0f\t\t%6.4f\t\t%6.4f\n',m,p(m),Np(m),C(m),E(m))
end
fprintf('\twAVG\t\t\t\t\t\t%6.4f\t\t%6.4f\n',sum(R(:,3).*C),sum(R(:,3).*E))

load('Source_Optimisation_Seizure');
fprintf('\n\n\t\t\t\t\tSeizure Results\n');
fprintf('\t  #\t\t\tp\t\tN_p\t\t  Corr\t\t Error\n');
for m = 1:M
    fprintf('\t  %1.0f\t\t\t%1.0f\t\t%1.0f\t\t%6.4f\t\t%6.4f\n',m,p(m),Np(m),C(m),E(m))
end
fprintf('\twAVG\t\t\t\t\t\t%6.4f\t\t%6.4f\n',sum(R(:,2).*C),sum(R(:,2).*E))

%% Load Optimal source positions
m = 1; % Patient 1 data
load(['Data\Optimisation\Positions\Optimal Positions\Patient_' num2str(m)]);
Ns = length(opt_pos_seiz);
ns = size(opt_pos_seiz{1},1);
Nb = length(opt_pos_back);
nb = size(opt_pos_back{1},1);
Xs = zeros(ns, Ns);
Ys = zeros(ns, Ns);
Zs = zeros(ns, Ns);
Xb = zeros(nb, Nb);
Yb = zeros(nb, Nb);
Zb = zeros(nb, Nb);
for ns = 1:Ns
    Xs(:,ns) = opt_pos_seiz{1,ns}(:,1);
    Ys(:,ns) = opt_pos_seiz{1,ns}(:,2);
    Zs(:,ns) = opt_pos_seiz{1,ns}(:,3);
end
[THI_s, PHI_s, R_s] = cart2sph(Xs,Ys,Zs);
for nb = 1:Nb
    Xb(:,nb) = opt_pos_back{1,nb}(:,1);
    Yb(:,nb) = opt_pos_back{1,nb}(:,2);
    Zb(:,nb) = opt_pos_back{1,nb}(:,3);
end
[THI_b, PHI_b, R_b] = cart2sph(Xb,Yb,Zb);

%% Plotting source posistions in 3D
p = sphere_cross(brain_R, pi/2, -pi/2, 32, 32);
figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(0,90); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xs,Ys,Zs,'.k','markersize',20); view(0,90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(90,0); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xs,Ys,Zs,'.k','markersize',20); view(90,0);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(0,0); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xs,Ys,Zs,'.k','markersize',20); view(0,0);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(0,90); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xb,Yb,Zb,'.b','markersize',20); view(0,90);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(90,0); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xb,Yb,Zb,'.b','markersize',20); view(90,0);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

figure('Color',[1,1,1]);
surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',0.5);
axis square; view(3); hold on; axis off;
lightangle(0,0); lighting phong; axis([-brain_R brain_R -brain_R brain_R -brain_R brain_R]);
plot3(Xb,Yb,Zb,'.b','markersize',20); view(0,0);
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'seizure_source_positions_xy','-dpdf','-opengl','-r400');
    print(2,'seizure_source_positions_yz','-dpdf','-opengl','-r400');
    print(3,'seizure_source_positions_xz','-dpdf','-opengl','-r400');
    print(4,'backrgound_source_positions_xy','-dpdf','-opengl','-r400');
    print(5,'backrgound_source_positions_yz','-dpdf','-opengl','-r400');
    print(6,'backrgound_source_positions_xz','-dpdf','-opengl','-r400');
end