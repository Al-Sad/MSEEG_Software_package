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
%                        Plotting MultiChannel Data
%
% Syntax : [out,tag] = plot_multichannel(x, space, fs, c, lw, legd, loc);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% x      : Multichannel Data (ch_n x N).
% space  : Vertical spacing between multichannel signals.
% fs     : Sampling frequency (Hz).
% c      : Line colour. see help of 'plot' (default : 'b').
% lw     : Line width. see help of 'plot' (default : 1).
% legd   : A flag to include the legend or not (default : 0).
% loc    : Legend location (default : 'NorthEastOutside').
%
% <OUTPUTs>
% out    : Plot handle.
% tag    : Channel labels.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE>
% x = randn(30,512);
% [h,tag] = plot_multichannel(x, 5, 1,'b');
% set(gca,'yticklabel',tag);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [out, tag] = plot_multichannel(x, space, fs, c, lw, legd, loc)
if(nargin < 4)
    c = 'b'; lw = 1; legd = 0;
elseif(nargin < 5)
    lw = 1; legd = 0;
elseif(nargin < 6)
    legd = 0;
elseif(nargin < 7 && legd == 1)
    loc = 'NorthEastOutside';
end
[ch_n, N] = size(x);
t = 0:1/fs:(N-1)/fs;
tag = cell(1,ch_n);
for i = 1:ch_n
    dat = x(i,:);
    out = plot(t, dat - space.*(i-1),'color',c,'linewidth',lw); hold on
    if(i < 10)
        tag{1,i} = ['CH. 0' num2str(i)];
    else
        tag{1,i} = ['CH. ' num2str(i)];
    end
end
set(gca,'yticklabel',{},'Ytick',(space*(1-ch_n)):space:0,'GridLineStyle','-');
grid on; ylim([(space*(1-ch_n))-space space]);
if(legd), legend(tag,'Location',loc);end
hold off
end