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
%                   Multi-sensor EEG propagation
%
% Syntax : [uni_Lamda,uni_Phi,bi_Lamda] = multi_sensor_propagation(event_p,fs,mP,v,varargin)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% event_p   : EEG source signals positions in 3D cartesian coordinates within
%             the brain sphere. It is a Px3 matrix, where P is the number of
%             EEG sources.
% fs        : Sampling frequency in Hz.
% mP        : Mean multi-sensor power of EEG.
% v         : Propagation speed of source signals in cm/s.
% varargin  : Optional input to hold custom sensor locations.
%
% <OUTPUTs>
% uni_Lamda : Multi-sensor unipolar attenuation matrix.
% uni_Phi   : Multi-sensor unipolar translation matrix.
% bi_Lamda  : Multi-sensor bipolar attenuation matrix.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% event_p = [-2.7 2.7 1.2 ; 0 -3 2 ; 4 0 0];
% fs = 32; mP = 100; v = 5;
% [uni_Lamda,uni_Phi,bi_Lamda] = multi_sensor_propagation(event_p,fs,mP,v);
% disp([uni_Lamda, uni_Phi]); disp(bi_Lamda);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [uni_Lamda,uni_Phi,bi_Lamda] = multi_sensor_propagation(event_p,fs,mP,v,varargin)
%% Parameters
[P, ~]   = size(event_p);       % Number of source signals
scalp_R  = 5.95;                % Head Radius in cm
th1      = 0.29;                % Scalp thickness in cm
th2      = 0.6;                 % Skull thickness in cm
th3      = 0.3;                 % CSF thickness in cm
% Spheres Radii (from center to the outer shell edge)
skull_R  = scalp_R - th1;       % Skull Radius in cm
csf_R    = skull_R - th2;       % CSF Radius in cm
brain_R  = csf_R - th3;         % Brain Radius in cm
% Tissue Optical Properties [brain, CSF, skull, scalp]
Rk  = [1.3,1.3,1.3,1.3];        % Refractive Indices
Ua  = [0.425,0.041, 0.16,0.18]; % Absorption Coefficients in 1/cm

%% EEG Sensors
if(isempty(varargin))
    % Unipolar EEG configuration
    [elec_p, untag] = electrode_10_20(scalp_R);
    [ch_n,~] = size(elec_p);        % Number of electrodes
    % Bipolar EEG configuration
    bitag = {'F_4','T_4';'T_4','T_6';'T_6','O_2';'F_3','T_3';'T_3','T_5';'T_5','O_1'
        'F_4','C_4';'C_4','P_4';'P_4','O_2';'F_3','C_3';'C_3','P_3';'P_3','O_1';'T_4','C_4'
        'C_4','C_z';'C_z','C_3';'C_3','T_3';'T_6','P_4';'P_4','P_z';'P_z','P_3';'P_3','T_5'};
    % Calculating Bipolar indicies
    [In, Im] = size(bitag);
    ind = zeros(In, Im);
    for m = 1:Im
        for n = 1:In
            for i = 1:length(untag)
                tmp = strcmp(bitag{n,m},untag{i});
                if(tmp)
                    ind(n,m) = i;
                end
            end
        end
    end
else
    elec_p = varargin{1};
    [ch_n,~] = size(elec_p);
end
%% Initialisation
lamda  = zeros(ch_n, P);
phi    = zeros(ch_n, P);

%% Main
for p = 1:P
    % Ghost Electrodes
    brain_ghost = line_sphere(event_p(p,:), elec_p, brain_R);
    csf_ghost   = line_sphere(event_p(p,:), elec_p, csf_R);
    skull_ghost = line_sphere(event_p(p,:), elec_p, skull_R);
    
    % Calculating signal path distances
    brain_dx = brain_ghost(:,1) - repmat(event_p(p,1), ch_n, 1);
    brain_dy = brain_ghost(:,2) - repmat(event_p(p,2), ch_n, 1);
    brain_dz = brain_ghost(:,3) - repmat(event_p(p,3), ch_n, 1);
    brain_D  = sqrt( (brain_dx).^2 + (brain_dy).^2 + (brain_dz).^2 ); % in cm
    
    csf_dx = csf_ghost(:,1) - brain_ghost(:,1);
    csf_dy = csf_ghost(:,2) - brain_ghost(:,2);
    csf_dz = csf_ghost(:,3) - brain_ghost(:,3);
    csf_D = sqrt( (csf_dx).^2 + (csf_dy).^2 + (csf_dz).^2 ); % in cm
    
    skull_dx = skull_ghost(:,1) - csf_ghost(:,1);
    skull_dy = skull_ghost(:,2) - csf_ghost(:,2);
    skull_dz = skull_ghost(:,3) - csf_ghost(:,3);
    skull_D = sqrt( (skull_dx).^2 + (skull_dy).^2 + (skull_dz).^2 ); % in cm
    
    scalp_dx = elec_p(:,1) - skull_ghost(:,1);
    scalp_dy = elec_p(:,2) - skull_ghost(:,2);
    scalp_dz = elec_p(:,3) - skull_ghost(:,3);
    scalp_D = sqrt( (scalp_dx).^2 + (scalp_dy).^2 + (scalp_dz).^2 ); % in cm
    
    D = scalp_D + skull_D + csf_D + brain_D;
    
    % Attenuation Factor Calculation
    Ao  = exp(Ua(2)*th3+Ua(3)*th2+Ua(4)*th1); % Initial Intensity
    R12 = ((Rk(1) - Rk(2))/(Rk(1) + Rk(2)))^2;
    R23 = ((Rk(2) - Rk(3))/(Rk(2) + Rk(3)))^2;
    R34 = ((Rk(3) - Rk(4))/(Rk(3) + Rk(4)))^2;
    RR  = (1-R12).*(1-R23).*(1-R34);
    lamda(:,p) = sqrt(Ao.*RR.*exp(-Ua(1)*brain_D-Ua(2)*csf_D-Ua(3)*skull_D-Ua(4)*scalp_D));
    
    % Phase Factor Calculation
    delay_temp = D./v;
    phi(:,p) = round(delay_temp*fs); % Quantifying the step with respect to fs
end
uni_Lamda = sum(lamda,2)./P;

%% Calculating bipolar attenuation and translation matrices
if P == 1
    uni_Phi = phi;
else
    uni_Phi = min(phi,[],2);
end
bi_Lamda = [];
if(isempty(varargin))
    bi_Lamda = zeros(In,1);
    for j = 1:In
        bi_Lamda(j,1) = uni_Lamda(ind(j,1),1) + uni_Lamda(ind(j,2),1);
    end
end
bi_Lamda = sqrt(mP).*(bi_Lamda./1.9263);

end