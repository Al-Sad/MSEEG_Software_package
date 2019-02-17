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
%                       Fractal signal generator
%
% Syntax : Fsig = FDsig(FD, fs, N);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% FD     : Fractal dimension. 
% fs     : Sampling frequency in Hz.
% N      : Number of samples.
%
% <OUTPUTs>
% Fsig   : Fractal signal.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% s = FDsig(2, 1, 256); plot(s);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Originally written by: Luke Rankine ~ October 2005.
% L. Rankine, N. Stevenson, M. Mesbah, B. Boashash, A Nonstationary
% Model of Newborn EEG, IEEE Transactions on Biomedical Engineering
% 54 (1) (2007) 19–28. doi:10.1109/TBME.2006.886667.
% Modified by: Mohammad Al-Sa'd ~ January 2019.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function Fsig = FDsig(FD, fs, N)
%% Parameters
gamma  = 5-2*FD;                   % Power Spectrum Power Law Exponent
Ncomp  = length(0:floor(N/2));     % Number of signal components
theta  = 2*pi*rand(1,N);           % Random phase factor
fp     = linspace(0,fs/2,Ncomp);   % Positive frequency vector
fn     = sort(-1*fp(2:Ncomp-1));   % Negative frequency vector
f      = [fp, fn];                 % Complete frequency vector

%% constructive and destructive interferences
sig = zeros(1,N);
for i = 2:Ncomp
    if i < Ncomp
        sig = sig + cos(2*pi*(i-1)/N *(0:N-1) + theta(i));
    else
        sig = sig + 0*cos(2*pi*(i-1)/N *(0:N-1)); % Removed Nyquist component
    end
end

%% Giving the Magnitude Spectrum the Correct Power Law
FT    = fft(sig);     % Fourier Transform of the accumulated signal
PS    = abs(FT).^2;   % Magnitude Power Spectrum
PSm   = [PS(1), zeros(1,length(PS)-1)];
for i = 2:length(PS)
   PSm(i) = (1/(abs(f(i))^gamma)) * PS(i);
end

%% Fourier Transform of the Fractal Signal
FT_new = zeros(1,length(FT));
for i = 1:length(FT)
    FT_new(i) = sqrt(PSm(i)) * exp(-1j*angle(FT(i)));
end

%% Time-Domain Fractal Signal
Fsig = real(ifft(FT_new));
end