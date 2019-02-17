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

%% Parameters
fs    = 32;        % Sampling frequency
R     = 5.95;      % Head Radius in cm
N     = 256;       % Number of samples for distribution
mP    = 100;       % Initial intensity
speed = 1;         % Propagation speed (cm/s)

%% Event positions as defined in Experiment 2 on pp 12
event_p = [-2.7 2.7 1.2 ; 0 -3 2 ; 4 0 0];

%% Main
for k = 1:length(event_p)
    % Power Dispersion
    [Vs, Phi, Scalp_p] = power_phase_dist(event_p(1:k,:), N, fs, mP, speed);
    % Plotting
    figure('Color',[1,1,1],'Position',[100 10 700 550]); colormap(jet);
    surface(Scalp_p{1}, Scalp_p{2}, Scalp_p{3}, Vs,'edgealpha',0); hold on
    [x,y] = draw_ellipse(0,R,1,0.5,N);
    plot3(x,y,zeros(N,1),'-k','linewidth',3); hold on
    [x,y] = draw_ellipse(0,-R,1,0.5,N);
    plot3(x,y,zeros(N,1),'-k','linewidth',3); hold on
    [x,y] = draw_ellipse(0,0,R,R,N);
    plot3(x,y,repmat(10,N,1),'-k','linewidth',3); hold on
    line([-R,-R-1], [-0.5, 0.01],'Color','k','linewidth',3); hold on;
    line([-R,-R-1], [0.5, -0.01],'Color','k','linewidth',3);
    axis square; axis([-R-1 R+1 -R-1 R+1]); view(2); axis off;
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
    
    figure('Color',[1,1,1],'Position',[100 10 700 550]); colormap(jet);
    surface(Scalp_p{1}, Scalp_p{2}, Scalp_p{3}, Phi,'edgealpha',0); hold on
    [x,y] = draw_ellipse(0,R,1,0.5,N);
    plot3(x,y,zeros(N,1),'-k','linewidth',3); hold on
    [x,y] = draw_ellipse(0,-R,1,0.5,N);
    plot3(x,y,zeros(N,1),'-k','linewidth',3); hold on
    [x,y] = draw_ellipse(0,0,R,R,N);
    plot3(x,y,repmat(10,N,1),'-k','linewidth',3); hold on
    line([-R,-R-1], [-0.5, 0.01],'Color','k','linewidth',3); hold on;
    line([-R,-R-1], [0.5, -0.01],'Color','k','linewidth',3);
    axis square; axis([-R-1 R+1 -R-1 R+1]); view(2); axis off;
    set(gcf,'Units','inches'); screenposition = get(gcf,'Position'); 
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
    % Saving
    opt = input('Do you want to save the results (Y/N)\n','s');
    if(opt == 'y' || opt == 'Y')
        print(1,['power_dispersion_' num2str(k)],'-dpdf','-r512');
        print(2,['phase_dispersion_' num2str(k)],'-dpdf','-r512');
        close all;
    elseif(opt == 'n' || opt == 'N')
        close all;
    else
        return
    end
end
