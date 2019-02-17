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
addpath(genpath('Data\3D'));

%% Parameters
N = 32; M = 32;    % Structure meshing size
scalp_R = 5.95;    % Head Radius in cm
th1 = 0.29;        % Scalp thickness in cm
th2 = 0.6;         % Skull thickness in cm
th3 = 0.3;         % CSF thickness in cm

%% Main
% 3D newborn head
[F, V] = baby(scalp_R);
% Electrode Positions
[elec_p, tag] = electrode_10_20(scalp_R);
p = sphere_cross(0.3, 0, 2*pi, N, M);
sph_p = cell(length(elec_p),3);
for i = 1:length(elec_p)
    sph_p{i,1} = p{1} + elec_p(i,1);
    sph_p{i,2} = p{2} + elec_p(i,2);
    sph_p{i,3} = p{3} + elec_p(i,3);
end
% Head Model
scalp = sphere_cross(scalp_R, 0, 2*pi, N, M);
[ear_r, ear_l] = ellipsoid_model(scalp_R, 0, 2*pi, N, M);
nose = cone_model(scalp_R, 0, 2*pi, pi/3, 2, N/2, M/2);

%% Plotting: Figures 1 and 2
%%% Real Neonatal Head
figure('Color',[1,1,1],'Position',[100 100 500 500]);
trisurf(F, V(:,1), V(:,2), V(:,3),'Linestyle','none','FaceColor',[1,0.75,0.65])
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
axis equal; view(3); axis off
lightangle(0,90); lighting phong
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));
%%% Neonatal Head Model
figure('Color',[1,1,1],'Position',[100 100 500 500]);
surf(scalp{1}, scalp{2}, scalp{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_r{1}, ear_r{2}, ear_r{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_l{1}, ear_l{2}, ear_l{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(nose{1}, nose{2}, nose{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1);
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
axis equal; view(3); axis off
lightangle(0,90); lighting phong
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting: Figure 3
figure('Color',[1,1,1],'position',[100, 60, 1000, 600]);
ha = tight_subplot(1,2,[0 -0.2],[0 0.2],[0.01 0.01]);
axes(ha(1));
%%% Plotting Scalp Sphere Cross Section
for R1 = (scalp_R-th1):0.01:scalp_R
    p = sphere_cross(R1, -pi/6, pi, N, M);
    if(R1 ~= scalp_R && R1 ~= scalp_R-th1)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0,'FaceAlpha',1); hold on
    else
        a = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.5,'FaceAlpha',1); hold on
    end
end
pn = cone_model(scalp_R, 0, pi, pi/3, 2, N/2, M/2);
pr = ellipsoid_model(scalp_R, 0, 2*pi, N, M);
surf(pn{1}, pn{2}, pn{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(pr{1}, pr{2}, pr{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
%%% Plotting Skull Sphere Cross Section
for R2 = (scalp_R-th1-th2):0.01:(scalp_R-th1-0.01)
    p = sphere_cross(R2, -pi/2, pi, N, M);
    if(R2 ~= scalp_R-th1-0.01 && R2 ~= scalp_R-th1-th2)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 0],'EdgeAlpha',0,'FaceAlpha',0.8); hold on
    else
        b = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 0],'EdgeAlpha',0.8,'FaceAlpha',0.8); hold on
    end 
end
%%% Plotting CSF Sphere Cross Section
for R3 = (scalp_R-th1-th2 - th3):0.01:(scalp_R-th1-th2-0.01)
    p = sphere_cross(R3, -2*pi/3, pi, N, M);
    if(R3 ~= scalp_R-th1-th2-0.01 && R3 ~= scalp_R-th1-th2 - th3)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 1],'EdgeAlpha',0,'FaceAlpha',0.8); hold on
    else
        c = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 1],'EdgeAlpha',0.8,'FaceAlpha',0.8); hold on
    end
end
%%% Plotting Brain Sphere
R4 = (scalp_R-th1-th2 - th3) - 0.01;
p = sphere_cross(R4, 0, 2*pi, N, M);
d = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',1); hold on
hold off; axis equal; axis off
axis([-scalp_R-2.1 scalp_R -scalp_R scalp_R+0.5 -scalp_R scalp_R]); view(3)
lightangle(0,90); lighting phong
set(gca,'fontweight','bold','fontsize',16);
ll = legend([a, b],'Scalp: Thickness = 2.9 mm','Skull:  Thickness = 6 mm');
set(ll,'Position',[0.175500003620982 0.0222222222222222 0.298999992758036 0.0933333308498066]);

axes(ha(2));
%%% Plotting Scalp Sphere Cross Section
for R1 = (scalp_R-th1):0.01:scalp_R
    p = sphere_cross(R1, -pi/6, pi, N, M);
    if(R1 ~= scalp_R && R1 ~= scalp_R-th1)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0,'FaceAlpha',1); hold on
    else
        a = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.5,'FaceAlpha',1); hold on
    end
end
pn = cone_model(scalp_R, 0, pi, pi/3, 2, N/2, M/2);
pr = ellipsoid_model(scalp_R, 0, 2*pi, N, M);
surf(pn{1}, pn{2}, pn{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(pr{1}, pr{2}, pr{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
%%% Plotting Skull Sphere Cross Section
for R2 = (scalp_R-th1-th2):0.01:(scalp_R-th1-0.01)
    p = sphere_cross(R2, -pi/2, pi, N, M);
    if(R2 ~= scalp_R-th1-0.01 && R2 ~= scalp_R-th1-th2)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 0],'EdgeAlpha',0,'FaceAlpha',0.8); hold on
    else
        b = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 0],'EdgeAlpha',0.8,'FaceAlpha',0.8); hold on
    end 
end
%%% Plotting CSF Sphere Cross Section
for R3 = (scalp_R-th1-th2 - th3):0.01:(scalp_R-th1-th2-0.01)
    p = sphere_cross(R3, -2*pi/3, pi, N, M);
    if(R3 ~= scalp_R-th1-th2-0.01 && R3 ~= scalp_R-th1-th2 - th3)
        surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 1],'EdgeAlpha',0,'FaceAlpha',0.8); hold on
    else
        c = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[0, 1, 1],'EdgeAlpha',0.8,'FaceAlpha',0.8); hold on
    end
end
%%% Plotting Brain Sphere
R4 = (scalp_R-th1-th2-th3) - 0.01;
p = sphere_cross(R4, 0, 2*pi, N, M);
d = surf(p{1}, p{2}, p{3},'Linestyle','-','FaceColor',[1, 0, 0],'EdgeAlpha',0.8,'FaceAlpha',1); hold on
hold off; axis equal; axis off
axis([-scalp_R-3 scalp_R+3 -scalp_R-3 scalp_R+3]); view(2);
lightangle(0,90); lighting phong
set(gca,'fontweight','bold','fontsize',16);
legend([c, d],'CSF:   Thickness = 3 mm','Brain:  Radius = 47.6 mm','location','South');
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Plotting: Figure 4
figure('Color',[1,1,1],'Position',[100 60 1000 600]);
ha = tight_subplot(1,2,[0 -0.1],[0 0],[0 0]);
axes(ha(1));
surf(scalp{1}, scalp{2}, scalp{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_r{1}, ear_r{2}, ear_r{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_l{1}, ear_l{2}, ear_l{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(nose{1}, nose{2}, nose{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
sh = 0.95;
for i = 1:length(elec_p)
    if i>=1 && i<4
    a = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','k','Linestyle','none'); hold on
    text(elec_p(i,1)-sh,elec_p(i,2)+sh,elec_p(i,3)+5,tag(i),'FontWeight','bold','FontSize',20,'color','k'); hold on
    end
    if i>=4 && i<12
    b = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','r','Linestyle','none'); hold on
    text(elec_p(i,1)-sh,elec_p(i,2)+sh,elec_p(i,3)+5,tag(i),'FontWeight','bold','FontSize',20,'color','r'); hold on
    end
    if i>=12 && i<20
    c = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','b','Linestyle','none'); hold on
    text(elec_p(i,1)-sh,elec_p(i,2)+sh,elec_p(i,3)+5,tag(i),'FontWeight','bold','FontSize',20,'color','b'); hold on
    end
    if i>=20 && i<=21
    d = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','g','Linestyle','none'); hold on
    text(elec_p(i,1)-sh,elec_p(i,2)+sh,elec_p(i,3)+5,tag(i),'FontWeight','bold','FontSize',20,'color','g');
    end
end
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{},'FontWeight','bold', 'FontSize',20)
axis equal; view(2); 
axis([-9 9 -9 9 -9 9]); axis off
lightangle(0,90); lighting phong
ll = legend([a, b, c, d],(sprintf('Midline\nElectrodes')),(sprintf('Left Hemisphere\nElectrodes'))...
    ,(sprintf('Right Hemisphere\nElectrodes')),...
    (sprintf('Front & Back\nElectrodes')),'Orientation','horizontal','Location','South');

set(ll,'Position',[0.11866668159763 0.0372222223712339 0.736999984234572 0.0933333308498065],...
    'Orientation','horizontal');

axes(ha(2));
surf(scalp{1}, scalp{2}, scalp{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_r{1}, ear_r{2}, ear_r{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(ear_l{1}, ear_l{2}, ear_l{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
surf(nose{1}, nose{2}, nose{3},'Linestyle','-','FaceColor',[1,0.75,0.65],'EdgeAlpha',0.2,'FaceAlpha',1); hold on
for i = 1:length(elec_p)
    if i>=1 && i<4
    a = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','k','Linestyle','none'); hold on
    end
    if i>=4 && i<12
    b = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','r','Linestyle','none'); hold on
    end
    if i>=12 && i<20
    c = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','b','Linestyle','none'); hold on
    end
    if i>=20 && i<=21
    d = surf(sph_p{i,1},sph_p{i,2},sph_p{i,3},'FaceColor','g','Linestyle','none'); hold on
    end
end
set(gca,'xticklabel',{},'yticklabel',{},'zticklabel',{})
axis equal; view(3); axis off
lightangle(0,90); lighting phong
set(gcf,'Units','inches'); screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],'PaperSize',screenposition(3:4));

%% Saving
opt = input('Do you want to save the results (Y/N)\n','s');
if(opt == 'y' || opt == 'Y')
    print(1,'real_head','-dpng','-r512');
    print(2,'sim_head','-dpng','-r512');
    print(3,'four_sphere', '-dpng', '-r512');
    print(4,'electrodes_placed', '-dpng', '-r512');
else
    return
end
