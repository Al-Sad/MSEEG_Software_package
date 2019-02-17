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
clear; close all; clc;
addpath(genpath('Data\EEG'));

%% Parameters
M = 7; % Total number of patients

%% Main
N = zeros(M,3);
for m = 1:M
    % Loading Real EEG
    load(['Patient_' num2str(m) '_Data_corrected.mat']);
    [N(m,1), ~, ~] = size(clean_EEG);
    [N(m,2), ~, ~] = size(EEG_S);
    [N(m,3), ~, ~] = size(EEG_B);
end
Nt = sum(N);
R  = zeros(M,3);
for m = 1:M
    R(m,:) = N(m,:)./Nt;
end
save('Data\EEG\Ratios','R');

%% Display Results
fprintf('\t N_e \t\t\t N_s \t\t\t N_b\n');
for m = 1:M
    fprintf('%4.0f (%4.2f%%)\t %4.0f (%4.2f%%) \t %4.0f (%4.2f%%)\n',...
        N(m,1),100*R(m,1),N(m,2),100*R(m,2),N(m,3),100*R(m,3))
end
disp(newline);
disp('           Total Number of Segments');
out = sprintf('\t  %4.0f \t\t %4.0f (%4.2f%%) \t %4.0f (%4.2f%%)',...
    Nt(1),Nt(2),100*Nt(2)./Nt(1),Nt(3),100*Nt(3)./Nt(1));
disp(out)
