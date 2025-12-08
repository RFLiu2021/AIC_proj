function [N ,peakIdx, peakV] = findSpikes_model(trace,threshold)
% THE FUNCTION DETECTS TIME STAMPS OF ACTION POTENTIALS OF NEURONAL DATA.
%
%   Input: 
%       'trace': data, the membrane voltage array of a neuron, 1-dim array.
%       'Theshold1': the value set for spike detection.
%   Output:
%       N - the number of spikes detected
%       out1 - time stamps of spike PEAKS.
%
%   Example:
%       [numberOfSpikes, spikePeakT] = findSpikes(data,-10)
%
%   Coded by:  Dr. Ruifeng Liu (Sun-Yatsen University) on 20210809


if nargin < 2
    threshold =0;
end

idx_tem = trace(:,1) > threshold;
% len_trace = length(trace);
% d = zeros(len_trace,1);
% d(idx_tem) = 1;
spikeBoarder = diff(idx_tem);
spike_L = find(spikeBoarder==1);
spike_R = find(spikeBoarder==-1);
if length(spike_L) ~= length(spike_R)
    error('Met mistake in counting spikes.');
end

nSpike = length(spike_L);
peak_idx = [];
peakV = [];
for i = 1:nSpike
    d_tem = trace(spike_L(i):spike_R(i));
    [peakV(i),peak_idx(i)] = max(d_tem);
    peak_idx(i) = spike_L(i)+peak_idx(i)-1;
end

if isempty( peak_idx )
    N=0;
    peakIdx = nan;
else
    N=length(peak_idx) ;
    peakIdx=peak_idx;
end


end


