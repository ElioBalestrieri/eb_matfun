function significance = my_TFCE(cond1, cond2, x_time)
%% performs threshold free cluster enhancememnt

% find dimensions
ntp = length(cond1);

dim_subj = find(size(cond1)~=ntp);

while dim_subj~=1

    cond1 = cond1';
    cond2 = cond2';
    dim_subj = find(size(cond1)~=ntp);

end
dim_tp = find(size(cond1)==ntp);
nsubj = size(cond1, 1);

fprintf('\nFound %i timepoints, assuming that n timepoints exceed N part', ntp)
fprintf('\nSubj dimension: %i\nTp dimension: %i\n', dim_subj, dim_tp)

%% compute tval series
[tvals, TFCE_t] = local_ttest_timeseries(cond1, cond2, ntp, x_time);

%% compute permutations
nperm = 100;
% reproduce 
rng(0)

bigdataset = cat(dim_subj, cond1, cond2);

perms_TFCE = nan(nperm, ntp);

for iPerm = 1:nperm
    
    idx_s = randsample(nsubj*2, nsubj*2);
    shuffled_cond1 = bigdataset(idx_s(1:nsubj), :);
    shuffled_cond2 = bigdataset(idx_s(nsubj+1:end), :);

    [~, swap_TFCE] = local_ttest_timeseries(shuffled_cond1,...
        shuffled_cond2, ntp, x_time);

    perms_TFCE(iPerm, :) = swap_TFCE;
    
    fprintf('\nDone permutation %i out of %i', iPerm, nperm)
    
end

%% compute percentiles

t_95perc = prctile(perms_TFCE, 95, 1);

figure; hold on;
plot(TFCE_t, 'LineWidth', 2)
plot(t_95perc, 'k', 'LineWidth', 2)

foo = 1;


end

%% ################## LOCAL FUNCTIONS ##########################

function [tvals, TFCE_vals] = local_ttest_timeseries(cond1, cond2, ntp, x_time)

tvals = nan(ntp,1);

for iTp = 1:ntp
    
    group1 = cond1(:,iTp);
    group2 = cond2(:,iTp);
    
    [~, ~, ~, stat] = ttest(group1, group2, 'Tail', 'right');

    tvals(iTp) = stat.tstat;
    
end

% constants
dh = .01;
E = .001;
H = 2;

TFCE_vals = compute_TFCE(tvals, x_time, dh, E, H);


end

%% sandbox
% map_x = tvals>0;
% cluster_conv = conv(map_x, ones(2, 1), 'same');
% 
% 
% fsample = 512;
% 
% for iSample = 1:numel(cluster_conv)
%     
%     this_val = cluster_conv(iSample);
%     
%     if this_val == 0 || this_val == 1
%         
%         cluster_dur(iSample) = 0;
%         
%     else
%         
%         cluster_dur(iSample) = cluster_dur(iSample-1)+1/fsample;
%         
%     end
%     
% end
%         
% 
% 
% 
% plot(test)
% 
% 
% func_int = xtime
% 
% 

