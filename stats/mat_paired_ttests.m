function [tvals, p, varargout] = mat_paired_ttests(mat1, mat2, cfg)
% the function computes ttests between 3d matrices.
% features on 1 and 2 dimensions, repetitions on 3rd dimension
% since it's paired ttest, matrices must have same size.
% 
% the function allows to perform cluster permutation test as well.
%
%% necessary in general
% cfg.tails = 'left', 'right', 'both'
%
%% necessary for permutation (fields + examples)
% cfg.alphat = .05;
% cfg.clustalplha_thresh = .05;
% cfg.npermute = 500;
%
%% wanna be fancy?
% cfg.plot = true
%
%% synopsis
% [tvals, p, clusterstat, signcluster_mask] = mat_paired_ttests(mat1, mat2, cfg)

% copy, correct, do whatever you want with this function. But at least
% thank me, or buy me a coffee when you meet me.
% started by eb 28-Jan-2020     

%% ttest part
[tvals, p] = local_ttest(mat1, mat2, cfg);

%% cluster permutation part

% cluster statistics
[clustermap, clusterstat] = local_determine_clusters(tvals, p, cfg);

% permutation
rng(1) % for reproducibility
clusterstat = local_permute(mat1, mat2, clusterstat, cfg);

% define logical mask for significant clusters
signcluster_mask = local_mask(clusterstat, clustermap, cfg);

%% plots (if wanted)
if isfield(cfg, 'plot')
    
    if cfg.plot
        
        bin_clustmap = clustermap;
        bin_clustmap(bin_clustmap~=0) = 1;
        
        figure; 
        h1 = imagesc(tvals); hold on;
        imcontour(bin_clustmap, 'k');
        c = colorbar; 
        c.Label.String = 't values';
        title(sprintf('repeated measure ttest (%i participants)', size(mat1, 3)));
        alphaval = double(signcluster_mask);
        alphaval(~signcluster_mask) = .5;
        set(h1, 'AlphaData', alphaval)
        
    end
    
end

% varargout
varargout{1} = clusterstat;
varargout{2} = signcluster_mask;


end


%% ####################### LOCAL FUNCTIONS

function [tvals, p] = local_ttest(mat1, mat2, cfg)

n1 = size(mat1, 3); n2 = size(mat2, 3);
if n1 ~= n2 
    error('n repetitions mismatch')
end
    
d = mat1-mat2;
mean_d = mean(d, 3);
sd_d = std(d, [], 3);
SE = sd_d/sqrt(n1);

tvals = mean_d./SE;

p = tcdf(tvals, n1);

switch lower(cfg.tails)
    
    case 'left'
        
        p = p;         %#ok<ASGSL>

    case 'right'
        
        p = 1-p;
        
    case 'both'
        
        p(p>.5) = 1-p(p>.5);
        p = 2*p;
        
    otherwise
        
        error('n tails not specified correctly')
        
end

end

function [clustermap, clusterstat] = local_determine_clusters(tvals, p, cfg)

% determine pvalues below alpha
lgcl_sign = p<cfg.alphat;

% determine clustermap
clustermap = bwlabel(lgcl_sign);

% determine clusterstatistics
vect_clust_lab = unique(clustermap)';

clusterstat = nan(numel(vect_clust_lab)-1,2);

for iClust = vect_clust_lab(2:end) % count from label n2, 0 is the null label

    swap_lgcl = clustermap==iClust;
    clusterstat(iClust, 2) = sum(tvals(swap_lgcl));
    clusterstat(iClust, 1) = iClust;
        
end

% sort clusters according to magnitude
clusterstat = sortrows(clusterstat, 2, 'descend');


end

function clusterstat =  local_permute(mat1, mat2, clusterstat, cfg)

% find number of repetitions
nrep = size(mat1, 3);

% cat matrices along 4th dimension
bigmat = cat(4, mat1, mat2);

% start permutations
clust_stat_dist = nan(cfg.nperm, 1);
[swap_mat1, swap_mat2] = deal(nan(size(mat1)));

% old_diff = mean(mat1, 3)- mean(mat2, 3);
% figure; imagesc(old_diff); colorbar;
 

for iPerm = 1:cfg.nperm
    
    % shuffle labels, but leave intact the comparison within participants
    for iPart = 1:nrep
       
        shuffled_idxs = randsample(2, 2);
        swap_mat1(:, :, iPart) = squeeze(bigmat(:,:,iPart,shuffled_idxs(1))); 
        swap_mat2(:, :, iPart) = squeeze(bigmat(:,:,iPart,shuffled_idxs(2)));
        
    end
    
%     new_diff = mean(swap_mat1, 3)- mean(swap_mat2, 3);
%     figure; imagesc(new_diff);
%     ylim([min(old_diff(:)), max(old_diff(:))])
%     colorbar;
    
    [swap_tvals, swap_p] = local_ttest(swap_mat1, swap_mat2, cfg);
    [~, swap_clusterstat] = local_determine_clusters(swap_tvals, swap_p, cfg);
    
    if ~isempty(swap_clusterstat)
        clust_stat_dist(iPerm) = swap_clusterstat(1, 2); % get the highest value
    else
        clust_stat_dist(iPerm) = 0;
    end
    
    if mod(iPerm, 10) == 0 
        fprintf('\n%i permutations done (over %i)', iPerm, cfg.nperm)
    end
        
%     waitforbuttonpress
    
end

% determine empirical cdf
[mat_(:,1), mat_(:,2)] = ecdf(clust_stat_dist);

% determine p values
nclusts = size(clusterstat, 1);

for iClusts = 1:nclusts
    
    p_cl = 1-mat_(find(mat_(:,2)<clusterstat(iClusts,2), 1, 'last'), 1);
    clusterstat(iClusts, 3) = p_cl;

end


end

function signcluster_mask = local_mask(clusterstat, clustermap, cfg)

who_is_significant = find(clusterstat(:,3)<cfg.clustalplha_thresh)';

signcluster_mask = false(size(clustermap));

for iCluster = who_is_significant
    
    clustcode = clusterstat(iCluster, 1);
    signcluster_mask = signcluster_mask | clustermap == clustcode;
    
end

end