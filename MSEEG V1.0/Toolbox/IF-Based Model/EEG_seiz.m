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
%                      Neonatal Seizure EEG Simulator
%
% Syntax : sig = EEG_seiz(N, fs)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% N      : Number of samples.
% fs     : Sampling frequency in Hz.
%
% <OUTPUTs>
% sig    : Simulated neonatal seizure (abnormal) EEG.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% s = EEG_seiz(256, 20); plot(s);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Originally written by: Nathan Stevenson.
% L. Rankine, N. Stevenson, M. Mesbah, B. Boashash, A Nonstationary
% Model of Newborn EEG, IEEE Transactions on Biomedical Engineering
% 54 (1) (2007) 19–28. doi:10.1109/TBME.2006.886667.
% Modified by: Mohammad Al-Sa'd ~ January 2019.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function sig = EEG_seiz(N, fs)
%% Parameters
K    = 5;   % maximum number of harmonics or components
M    = 3;   % number of pieces in IF law
% Harmonic Ratios (Gain Factors)
R(1) = 1;                                % Constant
R(2) = 1.0*betarnd(1.7,3.2,1,1) + 0.2;   % Range: 0.2-1.2
R(3) = 0.8*betarnd(1.5,4.1,1,1) + 0.2;   % Range: 0.2-1.0
R(4) = 0.4*betarnd(1.9,3.6,1,1) + 0.2;   % Range: 0.2-0.6
R(5) = 0.2*betarnd(1.4,1.2,1,1) + 0.2;   % Range: 0.2-0.4
P    = round(7*betarnd(1.8,3.0,1,K)+1);  % Number of Turning Points
V    = betarnd(3.9,8.0,K,10);            % Normalized Variation
s    = 0.15*betarnd(69.1,69.8,1,M)-0.075; % slopes of IF law (Zeta)
fst  = lognrnd(-0.17,0.55,1,1)+0.425;    % starting frequency
tp   = 1+round((N-1)*rand(1,3));         % location of IF law turning points
phs  = 2*pi*rand(1,K)-pi;                % initial phase of harmonics for seizure

%% Checking and completing location of IF law turning points
ref = find(tp==1);
tp(ref) = tp(ref)+2;
clear ref
ref = find(tp==N);
tp(ref) = tp(ref)-2;
tp = [1 tp N];
ref1 = find(diff(tp) == 0);  % check for time duplicates
tp((ref1)) = tp((ref1))+2;
tp = sort(tp);

%% main
% GENERATE HARMONIC COMPONENT AMPLITUDES
AM = am_gen(K, R, P, V, N, fs);
% GENERATE PIECEWISE IF LAW
[t,ife] = iflaw_gen(tp, s, fst, fs);
% GENERATE SEIZURE
sig = seiz_gen(ife, AM, phs, fs);
sig = sig./sqrt(sum(sig.^2)/N);

%% Auxiliary functions
    function [t,ife] = iflaw_gen(tp, s, fst, fs)
        %  ifl = iflaw_gen(tp,sl,fst,N)
        %  The function generates a piecewise linear IF law for the seizure
        %  generation algorithm. Specifically, this IF law is constructed of 3
        %  pieces.
        %
        %  Inputs: tp = turning points in samples (M+1x1 vector)
        %          s = slopes in Hz/sec (Mx1 vector)
        %          fst = start frequency in Hz (scalar)
        %
        %  Output: ife = 3 piece linear IF law (1xN vector)
        %          t = time in secs (1xN vector)
        %
        %  where M=3.
        %
        %  Nathan Stevenson
        %  Perinatal Research Centre
        %  January 2006
        %  n.stevenson@qut.edu.au
        N = tp(end);
        t = 0:1/fs:(N-1)/fs;
        C1 = fst;
        ife1 = s(1).*t+C1;                           % Piece 1
        C2 = ife1(tp(2))-s(2)*t(tp(2));              % Piece 2
        ife2 = s(2).*t+C2;
        C3 = ife2(tp(3))-s(3)*t(tp(3));              % Piece 3
        ife3 = s(3).*t+C3;
        % Apply rectangular functions
        ife = [ife1(tp(1):tp(2)), ife2(tp(2)+1:tp(3)), ife3(tp(3)+1:tp(end))];
    end

    function [AM] = am_gen(K,R,P,V,N,fs)
        %  AM = am_gen(K,R,V,P);
        %
        %  The function generates the amplitude modulation for the harmonic
        %  components.
        %
        %  Inputs: K  = is the harmonic or component number (default = 6)
        %          R  = is the harmonic ratio (1xK vector)
        %          P  = is the AM turning point number (1 x K vector) range is [1,5]
        %          V  = is the random amplitude values (K x Q matrix) Q is max P+2
        %               which is 7
        %          N  = discrete seizure length (default = 256)
        %          Fs = sampling frequency (default = 20Hz)
        %
        %  Output: AM = amplutide modulation function for each harmonic or
        %               component (KxN matrix)
        %          t  = time in secs (1xN vector)
        %
        %  Nathan Stevenson
        %  Perinatal Research Centre
        %  January 2006
        %  n.stevenson@qut.edu.au
        
        t = 0:1/fs:(N-1)/fs;
        
        for ii = 1:K
            test = 1;
            
            while test == 1
                v = R(ii)*(0.67+V(ii,1:P(1)+2));              % select amplitudes
                p=0:P(1)-1;
                q = 1+round(((N-1).*(p+rand(1,P(1)))./P(1))); % generate time indices
                ref = find(q==1);
                q(ref) = q(ref)+2;
                clear ref
                ref = find(q==N);
                q(ref) = q(ref)-2;
                q = [1 q N];
                
                ref1 = 1;
                while isempty(ref1)==0
                    [z, jj] = sort(q);
                    ref1 = min(find(diff(z) == 0));               % check for time duplicates
                    if q(jj(ref1))<=3
                        q(jj(ref1)) = q(jj(ref1))+2;
                    else
                        q(jj(ref1)) = q(jj(ref1))-2;
                    end
                    clear z jj
                end
                
                AM(ii,:) = spline(t(q),[0 v 0], t);           % perform interpolation
                % zero derivative ends
                AM(ii,(AM(ii,:)<0)) = 0;
                AM_test = max(AM(ii,:))/R(ii);
                if AM_test <= 1.67 && AM_test >= 0.67 
                    test = 0;
                else
                    test = 1;
                end
                
            end
            clear v p q r ref ref1 ref2 test AM_test
            
        end
    end

    function [sig] = seiz_gen(fi, AM, phs, fs)
        %  sig = seiz_gen(N, IF, AM, phs);
        %
        %  The function generates a seizure of length N using the if law IF and the
        %  amplitude modulation AM with initial phase phs.
        %
        %  Inputs: N  = is the seizure length
        %          IF = is the piecewise IF law
        %          AM  = is the amplitude modulation functions for the harmonics
        %          phs  = is the initial phase of each harmonic
        %
        %  Output: sig = newborn EEG seizure signal
        %
        %  Nathan Stevenson
        %  Perinatal Research Centre
        %  January 2006
        %  n.stevenson@qut.edu.au
        % Fundamental IF's and harmonics (normalised)
        sig = 0;
        for k = 1:K
            ffi = k.*fi;
            b   = AM(k,:);
            b(ffi>fs/2) = 0;
            sig = sig + b.*sin(2*pi.*cumtrapz(ffi)/fs+phs(k));
        end
    end
end