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
%                     Neonatal Background EEG Simulator
%
% Syntax : sig = EEG_back(N, fs)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N      : Number of samples.
% fs     : Sampling frequency in Hz.
%
% <OUTPUTs>
% sig    : Simulated neonatal background (Normal) EEG.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% s = EEG_back(256, 20); plot(s);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Notes : This simulator relies on betarnd to create the beta distribution.
%         Therefore, the statistics toolbox must be installed for this
%         simulator to function. This function may be edited by using a
%         gaussian distribution with the same mean and variance in place
%         of the beta distribution if the statistics toolbox is not
%         installed with similar results.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Originally written by: Luke Rankine ~ October 2005.
% L. Rankine, N. Stevenson, M. Mesbah, B. Boashash, A Nonstationary
% Model of Newborn EEG, IEEE Transactions on Biomedical Engineering
% 54 (1) (2007) 19–28. doi:10.1109/TBME.2006.886667.
% Modified by: Mohammad Al-Sa'd ~ January 2019.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function sig = EEG_back(N, fs)

%% Parameters
% The Fractal Dimension estimate is varied according to a Beta
% Distribution, with parameters, Alpha = 7.82 and Beta = 7.44
alpha   = 7.82; beta = 7.44;
FD      = 1 + betarnd(alpha, beta, 1, 1);
% Filter with random cutoffs
Wp_vect = 0.2*rand(1, 1) + 0.4;
Rs_vect = 7*rand(1, 1) + 6;

%% Main
Fsig = zeros(1,N);
% Simulating the constructive and destructive interference
for i = 1:15, Fsig = Fsig + FDsig(FD,fs,N); end

Wp = Wp_vect;
Ws = Wp -0.5*Wp;
Wp = Wp/(fs/2);
Ws = Ws/(fs/2);
Rp = 1;
Rs = Rs_vect;
[Nn, Wn] = buttord(Wp,Ws,Rp,Rs);
[Bb,Aa]  = butter(Nn,Wn,'high');
sig = filter(Bb,Aa,Fsig);
sig = sig./sqrt(sum(sig.^2)/N);
end