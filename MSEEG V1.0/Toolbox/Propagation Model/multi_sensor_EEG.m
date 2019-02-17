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
%                   Multi-sensor newborn EEG simulator
%
% Syntax : [X,Xb,Xs] = multi_sensor_EEG(mask,fs,event_b,event_s,mP_b,mP_s,v)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <INPUTs>
% mask    : A user-defined binary mask containing ones when seizure is
%           required and zeros when it is not.
% fs      : Sampling frequency in Hz.
% event_b : Background EEG source signals positions in 3D cartesian coordinates
%           within the brain sphere (default : event_b{1} = [0 0 0]).
%           "event_b" is a cell matrix that can contain a single stationary
%           position or multiple non-stationary positions, e.g.
%           "event_b{1} = [0 0 0]; event_b{2} = [0 1 0];"
%           Note that the number of background source positions must equal
%           the number of background segments in "mask", where a segment
%           contains 256 samples. 
% event_s : Seizure EEG source signals positions in 3D cartesian coordinates
%           within the brain sphere (default : event_s{1} = [3 0 0]).
%           "event_s" is a cell matrix that can contain a single stationary
%           position or multiple non-stationary positions, e.g.
%           "event_s{1} = [3 0 0]; event_s{2} = [0 3 0];"
%           Note that the number of seizure source positions must equal the
%           number of seizure segments in "mask", where a segment contains
%           256 samples. 
% mP_b    : Mean multi-sensor power of background EEG (default : mP_b = 166.3).
% mP_s    : Mean multi-sensor power of seizure EEG (default : mP_s = 2480.1).
% v       : Propagation speed of source signals in cm/s (default : v = 2).
%
% <OUTPUTs>
% X       : Complete simulated multi-sensor newborn EEG waveform containing
%           seizure and background patterns according to the defined binary
%           mask.
% Xb      : Simulated multi-sensor newborn background EEG waveform
%           generated according to the defined binary mask.
% Xs      : Simulated multi-sensor newborn seizure EEG waveform generated
%           according to the defined binary mask.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE #1>
% mask = [zeros(1,128) ones(1,256) zeros(1,128)];
% fs = 32;
% [X, Xb, Xs] = multi_sensor_EEG(mask, fs);
% figure; [~, tag] = plot_multichannel(Xb, 50, fs);
% title('Simulated multi-sensor newborn background EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(Xs, 50, fs);
% title('Simulated multi-sensor newborn seizure EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(X, 50, fs);
% title('Complete simulated multi-sensor newborn EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE #2>
% mask = [zeros(1,20) ones(1,300) zeros(1,250)];
% fs = 64;
% event_b{1} = event_pos(4);
% event_s{1} = event_pos(2);
% mP_b = 300;
% mP_s = 1000;
% [X, Xb, Xs] = multi_sensor_EEG(mask,fs,event_b,event_s,mP_b,mP_s);
% figure; [~, tag] = plot_multichannel(Xb, 50, fs);
% title('Simulated multi-sensor newborn background EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(Xs, 50, fs);
% title('Simulated multi-sensor newborn seizure EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(X, 50, fs);
% title('Complete simulated multi-sensor newborn EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% <EXAMPLE #3>
% mask = [zeros(1,256) ones(1,512) zeros(1,256)];
% fs = 32;
% event_b{1,1} = event_pos(3);
% event_b{1,2} = event_pos(3);
% event_s{1,1} = event_pos(1);
% event_s{1,2} = event_pos(3);
% mP_b = 300;
% mP_s = 1000;
% [X, Xb, Xs] = multi_sensor_EEG(mask,fs,event_b,event_s,mP_b,mP_s);
% figure; [~, tag] = plot_multichannel(Xb, 50, fs);
% title('Simulated multi-sensor newborn background EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(Xs, 50, fs);
% title('Simulated multi-sensor newborn seizure EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
% figure; [~, tag] = plot_multichannel(X, 50, fs);
% title('Complete simulated multi-sensor newborn EEG');
% xlabel('Time (s)'); set(gca,'yticklabel',tag);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

function [X, Xb, Xs] = multi_sensor_EEG(mask, fs, event_b, event_s, mP_b, mP_s, v)
%% Constants
scalp_R = 5.95;     % Head Radius in cm

%% Main input checking
if (nargin < 2 || ~exist('mask','var') || ~exist('fs','var') || isempty(mask) || isempty(fs))
    error('ERROR: This function requires atleast a mask and a sampling frequency ...');
end
if(length(mask) == 1)
    error('ERROR: Mask length should corrospond to the signal length. It has to be more than 1 ...');
end
if(length(fs) > 1 || fs <= 0)
    error('ERROR: Sampling frequency has to be a positive number ...')
end
if(sum(mask == 1 | mask == 0) ~= length(mask))
    error('ERROR: Mask should only contain ones and zeros ...')
end

%% Auxillary input checking
if(~exist('event_s','var') || isempty(event_s))
    event_s{1,1} = [3 0 0];
end
if(~exist('event_b','var') || isempty(event_b))
    event_b{1,1} = [0 0 0];
end
if(~iscell(event_b) || ~iscell(event_s))
    error('ERROR: Source positions must be a cell of 3D positions ...');
end
if(~exist('mP_b','var') || isempty(mP_b))
    mP_b = 166.3;
end
if(~exist('mP_s','var') || isempty(mP_s))
    mP_s = 2480.1;
end
if(~exist('speed','var') || isempty(v))
    v = 2;
end

%% EEG Sensors
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

%% Initialisation
N  = length(mask);
Xb = zeros(In, N);
Xs = zeros(In, N);
Ts = 8*fs;

%% Checking source positions
Ns = length(event_s);
Nb = length(event_b);
if(sum(mask) == N) % Mask is all ones
    if(Ns > 1 && ~(Ns == ceil(N/Ts)))
        error('ERROR: Seizure source positions must be equal to the segment number ...');
    end
    cnt1 = 0; cnt2 = 0; flag = 1;
    for i = 1:N
        if(mask(i) == 1 && flag == 1)
            cnt1 = cnt1 + 1;
            tis(1,cnt1) = i;
            flag = 0;
        elseif(mask(i) == 0 && flag == 0)
            cnt2 = cnt2 + 1;
            tfs(1,cnt2) = i-1;
            flag = 1;
        elseif(mask(i) == 1 && i == N)
            cnt2 = cnt2 + 1;
            tfs(1,cnt2) = N;
        end
    end
    ns = length(tis);
    Ls = [0 tfs-tis+1];
elseif(sum(mask) == 0) % Mask is all zeros
    if(Nb > 1 && ~(Nb == ceil(N/Ts)))
        error('ERROR: Background source positions must be equal to the segment number ...');
    end
    cnt1 = 0; cnt2 = 0; flag = 1;
    for i = 1:N
        if(mask(i) == 0 && flag == 1)
            cnt1 = cnt1 + 1;
            tib(1,cnt1) = i;
            flag = 0;
        elseif(mask(i) == 1 && flag == 0)
            cnt2 = cnt2 + 1;
            tfb(1,cnt2) = i-1;
            flag = 1;
        elseif(mask(i) == 0 && i == N)
            cnt2 = cnt2 + 1;
            tfb(1,cnt2) = N;
        end
    end
    nb = length(tib);
    Lb = [0 tfb-tib+1];
elseif(sum(mask)>0 && sum(mask)<N) % Mask is alternating
    cnt1 = 0; cnt2 = 0; flag = 1;
    for i = 1:N
        if(mask(i) == 1 && flag == 1)
            cnt1 = cnt1 + 1;
            tis(1,cnt1) = i;
            flag = 0;
        elseif(mask(i) == 0 && flag == 0)
            cnt2 = cnt2 + 1;
            tfs(1,cnt2) = i-1;
            flag = 1;
        elseif(mask(i) == 1 && i == N)
            cnt2 = cnt2 + 1;
            tfs(1,cnt2) = N;
        end
    end
    ns = length(tis);
    Ls = [0 tfs-tis+1];
    if(Ns > 1 && ~(Ns == sum(ceil(Ls/Ts))))
        error('ERROR: Seizure source positions must be equal to the segment number ...');
    end
    cnt1 = 0; cnt2 = 0; flag = 1;
    for i = 1:N
        if(mask(i) == 0 && flag == 1)
            cnt1 = cnt1 + 1;
            tib(1,cnt1) = i;
            flag = 0;
        elseif(mask(i) == 1 && flag == 0)
            cnt2 = cnt2 + 1;
            tfb(1,cnt2) = i-1;
            flag = 1;
        elseif(mask(i) == 0 && i == N)
            cnt2 = cnt2 + 1;
            tfb(1,cnt2) = N;
        end
    end
    nb = length(tib);
    Lb = [0 tfb-tib+1];
    if(Nb > 1 && ~(Nb == sum(ceil(Lb/Ts))))
        error('ERROR: Background source positions must be equal to the segment number ...');
    end
end

%% Checking if positions are stationary non-stationary
Ps  = zeros(In,Ns);
Phi = zeros(ch_n,Ns);
Pb  = zeros(In,Nb);
if(Ns == 1)
    [~, Phi, Ps] = multi_sensor_propagation(event_s{1,1}, fs, mP_s, v);
else
    for i = 1:Ns
        [~, Phi(:,i), Ps(:,i)] = multi_sensor_propagation(event_s{1,i}, fs, mP_s, v);
    end
end
if(Nb == 1)
    [~, ~, Pb] = multi_sensor_propagation(event_b{1,1}, fs, mP_b, v);
else
    for i = 1:Nb
        [~, ~, Pb(:,i)] = multi_sensor_propagation(event_b{1,i}, fs, mP_b, v);
    end
end

%% Generate multi-sensor EEG
if(sum(mask) == N) % Mask is all ones
    dump = EEG_seiz(N, fs);
    % Overall template
    Min_phi = min(Phi(:,1));
    Max_phi = max(Phi(:,end));
    uni_phi_st = 1 + Phi(:,1)   - Min_phi;
    uni_phi_fi = N + Phi(:,end) - Max_phi;
    for j = 1:In
        bi_phi_st = min(uni_phi_st(ind(j,1)),uni_phi_st(ind(j,2)));
        bi_phi_fi = max(uni_phi_fi(ind(j,1)),uni_phi_fi(ind(j,2)));
        Xs(j,bi_phi_st:bi_phi_fi) = dump(1,1:bi_phi_fi-bi_phi_st+1);
    end
    % Single sources in time
    if(Ns == 1)
        % Integer number of segments
        if(~rem(N,Ts))
            % Adjust powers according to positions
            for i = 1:N/Ts
                EEG_st = 1  + (i-1)*Ts;
                EEG_fi = Ts + (i-1)*Ts;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                    end
                end
            end
            % Floating number of segments
        else
            remaining = rem(N,Ts);
            K = N - remaining;
            % Adjust powers according to positions
            for i = 1:K/Ts
                EEG_st = 1  + (i-1)*Ts;
                EEG_fi = Ts + (i-1)*Ts;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                    end
                end
            end
            EEG_st = 1 + EEG_fi;
            EEG_fi = N;
            for j = 1:In
                if(sum(Xs(j,EEG_st:EEG_fi).^2))
                    Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                end
            end
        end
        % Multi-sources in time
    else
        % Integer number of segments
        if(~rem(N,Ts))
            % Adjust powers according to positions
            for i = 1:N/Ts
                EEG_st = 1  + (i-1)*Ts;
                EEG_fi = Ts + (i-1)*Ts;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,i).*Xs(j,EEG_st:EEG_fi);
                    end
                end
            end
            % Floating number of segments
        else
            remaining = rem(N,Ts);
            K = N - remaining;
            % Adjust powers according to positions
            for i = 1:K/Ts
                EEG_st = 1  + (i-1)*Ts;
                EEG_fi = Ts + (i-1)*Ts;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,i).*Xs(j,EEG_st:EEG_fi);
                    end
                end
            end
            EEG_st = 1 + EEG_fi;
            EEG_fi = N;
            for j = 1:In
                if(sum(Xs(j,EEG_st:EEG_fi).^2))
                    Xs(j,EEG_st:EEG_fi) = Ps(j,end).*Xs(j,EEG_st:EEG_fi);
                end
            end
        end
    end
    % Add background to delayed segments
    bi_phi_st = min(uni_phi_st(ind(:,1)),uni_phi_st(ind(:,2)));
    bi_phi_fi = max(uni_phi_fi(ind(:,1)),uni_phi_fi(ind(:,2)));
    Sb = sum(bi_phi_st-1);
    Fb = sum(N-bi_phi_fi);
    temp1 = EEG_back(Sb +(rem(Sb,2)), fs);
    temp2 = EEG_back(Fb +(rem(Fb,2)), fs);
    for j = 1:In
        st_s = 1 + sum(bi_phi_st(1:(j-1))-1);
        fi_s = sum(bi_phi_st(1:j)-1);
        st_f = 1 + sum(N-bi_phi_fi(1:(j-1)));
        fi_f = sum(N-bi_phi_fi(1:j));
        Xb(j,1:(bi_phi_st(j)-1)) = Pb(j,1).*temp1(1,st_s:fi_s);
        Xb(j,(bi_phi_fi(j)+1):N) = Pb(j,1).*temp2(1,st_f:fi_f);
    end
    X = Xb + Xs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(sum(mask) == 0) % Mask is all zeros
    % Single sources in time
    if(Nb == 1)
        % Integer number of segments
        if(~rem(N,Ts))
            % Adjust powers according to positions
            for i = 1:N/Ts
                temp = zeros(In,Ts);
                for j = 1:In
                    temp(j,:) = Pb(j,1).*EEG_back(Ts, fs);
                end
                st = 1  + (i-1)*Ts;
                fi = Ts + (i-1)*Ts;
                Xb(:,st:fi) = temp;
            end
            % Floating number of segments
        else
            remaining = rem(N,Ts);
            K = N - remaining;
            % Adjust powers according to positions
            for i = 1:K/Ts
                temp = zeros(In,Ts);
                for j = 1:In
                    temp(j,:) = Pb(j,1).*EEG_back(Ts, fs);
                end
                st = 1  + (i-1)*Ts;
                fi = Ts + (i-1)*Ts;
                Xb(:,st:fi) = temp;
            end
            temp = zeros(In,remaining);
            if(rem(remaining,2))
                for j = 1:In
                    tmp = Pb(j,end).*EEG_back(remaining+1, fs);
                    temp(j,:) = tmp(1,1:end-1);
                end
            else
                for j = 1:In
                    temp(j,:) = Pb(j,end).*EEG_back(remaining, fs);
                end
            end
            Xb(:,K+1:N) = temp;
        end
        
        % Multi-sources in time
    else
        % Integer number of segments
        if(~rem(N,Ts))
            % Adjust powers according to positions
            for i = 1:N/Ts
                temp = zeros(In,Ts);
                for j = 1:In
                    temp(j,:) = Pb(j,i).*EEG_back(Ts, fs);
                end
                st = 1 + (i-1)*Ts;
                fi = Ts + (i-1)*Ts;
                Xb(:,st:fi) = temp;
            end
            % Floating number of segments
        else
            remaining = rem(N,Ts);
            K = N - remaining;
            % Adjust powers according to positions
            for i = 1:K/Ts
                temp = zeros(In,Ts);
                for j = 1:In
                    temp(j,:) = Pb(j,i).*EEG_back(Ts, fs);
                end
                st = 1 + (i-1)*Ts;
                fi = Ts + (i-1)*Ts;
                Xb(:,st:fi) = temp;
            end
            temp = zeros(In,remaining);
            if(rem(remaining,2))
                for j = 1:In
                    tmp = Pb(j,end).*EEG_back(remaining+1, fs);
                    temp(j,:) = tmp(1,1:end-1);
                end
            else
                for j = 1:In
                    temp(j,:) = Pb(j,end).*EEG_back(remaining, fs);
                end
            end
            Xb(:,K+1:N) = temp;
        end
    end
    X = Xb + Xs;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif(sum(mask)>0 && sum(mask)<N) % Mask is alternating
    SS = Ls - rem(Ls,Ts) + Ts;
    dump = EEG_seiz(sum(SS(2:end)), fs);
    % Single sources in time
    if(Ns == 1)
        LL = Ls; LL(1) = Ts;
        for nns = 1:ns
            % Integer number of segments
            if(~rem(Ls(nns+1),Ts))
                % Overall template
                Min_phi = min(Phi(:,1));
                Max_phi = max(Phi(:,end));
                uni_phi_st = tis(nns) + Phi(:,1)   - Min_phi;
                uni_phi_fi = tfs(nns) + Phi(:,end) - Max_phi;
                dump_uni_phi_st = 1 + sum(Ls(1:nns)) + Phi(:,1) - Min_phi;
                dump_uni_phi_fi = sum(Ls(1:nns+1)) + Phi(:,end) - Max_phi;
                for j = 1:In
                    bi_phi_st = min(uni_phi_st(ind(j,1)),uni_phi_st(ind(j,2)));
                    bi_phi_fi = max(uni_phi_fi(ind(j,1)),uni_phi_fi(ind(j,2)));
                    dump_bi_phi_st = min(dump_uni_phi_st(ind(j,1)),dump_uni_phi_st(ind(j,2)));
                    dump_bi_phi_fi = max(dump_uni_phi_fi(ind(j,1)),dump_uni_phi_fi(ind(j,2)));
                    Xs(j,bi_phi_st:bi_phi_fi) = dump(1,min(dump_uni_phi_st):...
                        (dump_bi_phi_fi-dump_bi_phi_st+min(dump_uni_phi_st)));
                end
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nns))/Ts:sum(Ls(1:nns+1))/Ts
                    cnt = cnt + 1;
                    EEG_st = tis(nns) + (cnt-1)*Ts;
                    EEG_fi = tis(nns) + cnt*Ts - 1;
                    for j = 1:In
                        if(sum(Xs(j,EEG_st:EEG_fi).^2))
                            Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                        end
                    end
                end
                % Add background to delayed segments
                bi_phi_st = min(uni_phi_st(ind(:,1)),uni_phi_st(ind(:,2)));
                bi_phi_fi = max(uni_phi_fi(ind(:,1)),uni_phi_fi(ind(:,2)));
                Sb = sum(bi_phi_st-tis(nns));
                Fb = sum(tfs(nns)-bi_phi_fi);
                temp1 = EEG_back(Sb +(rem(Sb,2)), fs);
                temp2 = EEG_back(Fb +(rem(Fb,2)), fs);
                for j = 1:In
                    st_s = 1 + sum(bi_phi_st(1:(j-1))-tis(nns));
                    fi_s = sum(bi_phi_st(1:j)-tis(nns));
                    st_f = 1 + sum(tfs(nns)-bi_phi_fi(1:(j-1)));
                    fi_f = sum(tfs(nns)-bi_phi_fi(1:j));
                    Xb(j,tis(nns):(bi_phi_st(j)-1)) = Pb(j,1).*temp1(1,st_s:fi_s);
                    Xb(j,(bi_phi_fi(j)+1):tfs(nns)) = Pb(j,1).*temp2(1,st_f:fi_f);
                end
                % Floating number of segments
            else
                remaining = rem(Ls(nns+1),Ts);
                Ls(nns+1) = Ls(nns+1) - remaining + Ts;
                LL = Ls; LL(1) = Ts;
                % Overall template
                Min_phi = min(Phi(:,1));
                Max_phi = max(Phi(:,end));
                uni_phi_st = tis(nns) + Phi(:,1)  - Min_phi;
                uni_phi_fi = tfs(nns) + Phi(:,end) - Max_phi;
                dump_uni_phi_st = 1 + sum(Ls(1:nns)) + Phi(:,1) - Min_phi;
                dump_uni_phi_fi = sum(Ls(1:nns+1)) + Phi(:,end) - Max_phi  + remaining - Ts;
                for j = 1:In
                    bi_phi_st = min(uni_phi_st(ind(j,1)),uni_phi_st(ind(j,2)));
                    bi_phi_fi = max(uni_phi_fi(ind(j,1)),uni_phi_fi(ind(j,2)));
                    dump_bi_phi_st = min(dump_uni_phi_st(ind(j,1)),dump_uni_phi_st(ind(j,2)));
                    dump_bi_phi_fi = max(dump_uni_phi_fi(ind(j,1)),dump_uni_phi_fi(ind(j,2)));
                    Xs(j,bi_phi_st:bi_phi_fi) = dump(1,min(dump_uni_phi_st):dump_bi_phi_fi-...
                        dump_bi_phi_st+min(dump_uni_phi_st));
                end
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nns))/Ts:sum(Ls(1:nns+1))/Ts-1
                    cnt = cnt + 1;
                    EEG_st = tis(nns) + (cnt-1)*Ts;
                    EEG_fi = tis(nns) + cnt*Ts - 1;
                    for j = 1:In
                        if(sum(Xs(j,EEG_st:EEG_fi).^2))
                            Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                        end
                    end
                end
                EEG_st = tis(nns) + (cnt)*Ts;
                EEG_fi = tis(nns) + (cnt)*Ts + remaining -1;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,1).*Xs(j,EEG_st:EEG_fi);
                    end
                end
                % Add background to delayed segments
                bi_phi_st = min(uni_phi_st(ind(:,1)),uni_phi_st(ind(:,2)));
                bi_phi_fi = max(uni_phi_fi(ind(:,1)),uni_phi_fi(ind(:,2)));
                Sb = sum(bi_phi_st-tis(nns));
                Fb = sum(tfs(nns)-bi_phi_fi);
                temp1 = EEG_back(Sb +(rem(Sb,2)), fs);
                temp2 = EEG_back(Fb +(rem(Fb,2)), fs);
                for j = 1:In
                    st_s = 1 + sum(bi_phi_st(1:(j-1))-tis(nns));
                    fi_s = sum(bi_phi_st(1:j)-tis(nns));
                    st_f = 1 + sum(tfs(nns)-bi_phi_fi(1:(j-1)));
                    fi_f = sum(tfs(nns)-bi_phi_fi(1:j));
                    Xb(j,tis(nns):(bi_phi_st(j)-1)) = Pb(j,1).*temp1(1,st_s:fi_s);
                    Xb(j,(bi_phi_fi(j)+1):tfs(nns)) = Pb(j,1).*temp2(1,st_f:fi_f);
                end
            end
        end
        % Multi-sources in time
    else
        Ls = [0 tfs-tis+1];
        LL = Ls; LL(1) = Ts;
        for nns = 1:ns
            % Integer number of segments
            if(~rem(Ls(nns+1),Ts))
                % Overall template
                Min_phi = min(Phi(:,sum(LL(1:nns))/Ts));
                Max_phi = max(Phi(:,sum(Ls(1:nns+1))/Ts));
                uni_phi_st = tis(nns) + Phi(:,sum(LL(1:nns))/Ts)  - Min_phi;
                uni_phi_fi = tfs(nns) + Phi(:,sum(Ls(1:nns+1))/Ts) - Max_phi;
                dump_uni_phi_st = 1 + sum(Ls(1:nns)) + Phi(:,sum(LL(1:nns))/Ts) - Min_phi;
                dump_uni_phi_fi = sum(Ls(1:nns+1)) + Phi(:,sum(Ls(1:nns+1))/Ts) - Max_phi;
                for j = 1:In
                    bi_phi_st = min(uni_phi_st(ind(j,1)),uni_phi_st(ind(j,2)));
                    bi_phi_fi = max(uni_phi_fi(ind(j,1)),uni_phi_fi(ind(j,2)));
                    dump_bi_phi_st = min(dump_uni_phi_st(ind(j,1)),dump_uni_phi_st(ind(j,2)));
                    dump_bi_phi_fi = max(dump_uni_phi_fi(ind(j,1)),dump_uni_phi_fi(ind(j,2)));
                    Xs(j,bi_phi_st:bi_phi_fi) = dump(1,min(dump_uni_phi_st):dump_bi_phi_fi-...
                        dump_bi_phi_st+min(dump_uni_phi_st));
                end
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nns))/Ts:sum(Ls(1:nns+1))/Ts
                    cnt = cnt + 1;
                    EEG_st = tis(nns) + (cnt-1)*Ts;
                    EEG_fi = tis(nns) + cnt*Ts - 1;
                    for j = 1:In
                        if(sum(Xs(j,EEG_st:EEG_fi).^2))
                            Xs(j,EEG_st:EEG_fi) = Ps(j,i).*Xs(j,EEG_st:EEG_fi);
                        end
                    end
                end
                % Add background to delayed segments
                bi_phi_st = min(uni_phi_st(ind(:,1)),uni_phi_st(ind(:,2)));
                bi_phi_fi = max(uni_phi_fi(ind(:,1)),uni_phi_fi(ind(:,2)));
                Sb = sum(bi_phi_st-tis(nns));
                Fb = sum(tfs(nns)-bi_phi_fi);
                temp1 = EEG_back(Sb +(rem(Sb,2)), fs);
                temp2 = EEG_back(Fb +(rem(Fb,2)), fs);
                for j = 1:In
                    st_s = 1 + sum(bi_phi_st(1:(j-1))-tis(nns));
                    fi_s = sum(bi_phi_st(1:j)-tis(nns));
                    st_f = 1 + sum(tfs(nns)-bi_phi_fi(1:(j-1)));
                    fi_f = sum(tfs(nns)-bi_phi_fi(1:j));
                    Xb(j,tis(nns):(bi_phi_st(j)-1)) = Pb(j,1).*temp1(1,st_s:fi_s);
                    Xb(j,(bi_phi_fi(j)+1):tfs(nns)) = Pb(j,1).*temp2(1,st_f:fi_f);
                end
                % Floating number of segments
            else
                remaining = rem(Ls(nns+1),Ts);
                Ls(nns+1) = Ls(nns+1) - remaining + Ts;
                LL = Ls; LL(1) = Ts;
                % Overall template
                Min_phi = min(Phi(:,sum(LL(1:nns))/Ts));
                Max_phi = max(Phi(:,sum(Ls(1:nns+1))/Ts));
                uni_phi_st = tis(nns) + Phi(:,sum(LL(1:nns))/Ts)  - Min_phi;
                uni_phi_fi = tfs(nns) + Phi(:,sum(Ls(1:nns+1))/Ts) - Max_phi;
                dump_uni_phi_st = 1 + sum(Ls(1:nns)) + Phi(:,sum(LL(1:nns))/Ts) - Min_phi;
                dump_uni_phi_fi = sum(Ls(1:nns+1)) + Phi(:,sum(Ls(1:nns+1))/Ts) - Max_phi  + remaining - Ts;
                for j = 1:In
                    bi_phi_st = min(uni_phi_st(ind(j,1)),uni_phi_st(ind(j,2)));
                    bi_phi_fi = max(uni_phi_fi(ind(j,1)),uni_phi_fi(ind(j,2)));
                    dump_bi_phi_st = min(dump_uni_phi_st(ind(j,1)),dump_uni_phi_st(ind(j,2)));
                    dump_bi_phi_fi = max(dump_uni_phi_fi(ind(j,1)),dump_uni_phi_fi(ind(j,2)));
                    Xs(j,bi_phi_st:bi_phi_fi) = dump(1,min(dump_uni_phi_st):dump_bi_phi_fi-...
                        dump_bi_phi_st+min(dump_uni_phi_st));
                end
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nns))/Ts:sum(Ls(1:nns+1))/Ts-1
                    cnt = cnt + 1;
                    EEG_st = tis(nns) + (cnt-1)*Ts;
                    EEG_fi = tis(nns) + cnt*Ts - 1;
                    for j = 1:In
                        if(sum(Xs(j,EEG_st:EEG_fi).^2))
                            Xs(j,EEG_st:EEG_fi) = Ps(j,i).*Xs(j,EEG_st:EEG_fi);
                        end
                    end
                end
                EEG_st = tis(nns) + (cnt)*Ts;
                EEG_fi = tis(nns) + (cnt)*Ts + remaining -1;
                for j = 1:In
                    if(sum(Xs(j,EEG_st:EEG_fi).^2))
                        Xs(j,EEG_st:EEG_fi) = Ps(j,i).*Xs(j,EEG_st:EEG_fi);
                    end
                end
                % Add background to delayed segments
                bi_phi_st = min(uni_phi_st(ind(:,1)),uni_phi_st(ind(:,2)));
                bi_phi_fi = max(uni_phi_fi(ind(:,1)),uni_phi_fi(ind(:,2)));
                Sb = sum(bi_phi_st-tis(nns));
                Fb = sum(tfs(nns)-bi_phi_fi);
                temp1 = EEG_back(Sb +(rem(Sb,2)), fs);
                temp2 = EEG_back(Fb +(rem(Fb,2)), fs);
                for j = 1:In
                    st_s = 1 + sum(bi_phi_st(1:(j-1))-tis(nns));
                    fi_s = sum(bi_phi_st(1:j)-tis(nns));
                    st_f = 1 + sum(tfs(nns)-bi_phi_fi(1:(j-1)));
                    fi_f = sum(tfs(nns)-bi_phi_fi(1:j));
                    Xb(j,tis(nns):(bi_phi_st(j)-1)) = Pb(j,1).*temp1(1,st_s:fi_s);
                    Xb(j,(bi_phi_fi(j)+1):tfs(nns)) = Pb(j,1).*temp2(1,st_f:fi_f);
                end
            end
        end
    end
    % Single sources in time
    if(Nb == 1)
        LL = Lb; LL(1) = Ts;
        for nnb = 1:nb
            % Integer number of segments
            if(~rem(Lb(nnb+1),Ts))
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nnb))/Ts:sum(Lb(1:nnb+1))/Ts
                    temp = zeros(In,Ts);
                    for j = 1:In
                        temp(j,:) = Pb(j,1).*EEG_back(Ts, fs);
                    end
                    cnt = cnt + 1;
                    EEG_st = tib(nnb) + (cnt-1)*Ts;
                    EEG_fi = tib(nnb) + cnt*Ts - 1;
                    Xb(:,EEG_st:EEG_fi) = temp;
                end
                % Floating number of segments
            else
                remaining = rem(Lb(nnb+1),Ts);
                Lb(nnb+1) = Lb(nnb+1) - remaining + Ts;
                LL = Lb; LL(1) = Ts;
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nnb))/Ts:sum(Lb(1:nnb+1))/Ts-1
                    temp = zeros(In,Ts);
                    for j = 1:In
                        temp(j,:) = Pb(j,1).*EEG_back(Ts, fs);
                    end
                    cnt = cnt + 1;
                    EEG_st = tib(nnb) + (cnt-1)*Ts;
                    EEG_fi = tib(nnb) + cnt*Ts - 1;
                    Xb(:,EEG_st:EEG_fi) = temp;
                end
                temp = zeros(In,remaining);
                for j = 1:In
                    temp(j,:) = Pb(j,end).*EEG_back(remaining, fs);
                end
                EEG_st = tib(nnb) + (cnt)*Ts;
                EEG_fi = tib(nnb) + (cnt)*Ts + remaining -1;
                Xb(:,EEG_st:EEG_fi) = temp;
            end
        end
        % Multi-sources in time
    else
        LL = Lb; LL(1) = Ts;
        for nnb = 1:nb
            % Integer number of segments
            if(~rem(Lb(nnb+1),Ts))
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nnb))/Ts:sum(Lb(1:nnb+1))/Ts
                    temp = zeros(In,Ts);
                    for j = 1:In
                        temp(j,:) = Pb(j,i).*EEG_back(Ts, fs);
                    end
                    cnt = cnt + 1;
                    EEG_st = tib(nnb) + (cnt-1)*Ts;
                    EEG_fi = tib(nnb) + cnt*Ts - 1;
                    Xb(:,EEG_st:EEG_fi) = temp;
                end
                % Floating number of segments
            else
                remaining = rem(Lb(nnb+1),Ts);
                Lb(nnb+1) = Lb(nnb+1) - remaining + Ts;
                LL = Lb; LL(1) = Ts;
                % Adjust powers according to positions
                cnt = 0;
                for i = sum(LL(1:nnb))/Ts:sum(Lb(1:nnb+1))/Ts-1
                    temp = zeros(In,Ts);
                    for j = 1:In
                        temp(j,:) = Pb(j,i).*EEG_back(Ts, fs);
                    end
                    cnt = cnt + 1;
                    EEG_st = tib(nnb) + (cnt-1)*Ts;
                    EEG_fi = tib(nnb) + cnt*Ts - 1;
                    Xb(:,EEG_st:EEG_fi) = temp;
                end
                temp = zeros(In,remaining);
                for j = 1:In
                    temp(j,:) = Pb(j,i).*EEG_back(remaining, fs);
                end
                EEG_st = tib(nnb) + (cnt)*Ts;
                EEG_fi = tib(nnb) + (cnt)*Ts + remaining -1;
                Xb(:,EEG_st:EEG_fi) = temp;
            end
        end
    end
    X = Xb + Xs;
end
end