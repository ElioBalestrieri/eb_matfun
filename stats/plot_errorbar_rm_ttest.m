function plot_errorbar_rm_ttest(matrix, labels)

hold on;

avg_ = mean(matrix);
err_ = std(matrix)/sqrt(size(matrix,1));

for iVect = 1:size(matrix, 2)
    this_vect = matrix(:,iVect);
    xsc_ = (randn(numel(this_vect),1))/30 + iVect;
    scatter(xsc_, this_vect, 20, [100 192 238]/255, 'filled')    
end

[~, p, ~, stat] = ttest(matrix(:, 1), matrix(:, 2));

stringsign1 = sprintf('t = %4.4f', stat.tstat);
stringsign2 = sprintf('p = %4.4f', p);

text(1.5-.15, mean(matrix(:)), {stringsign1, stringsign2})

errorbar(1:length(avg_), avg_, err_, '.k', 'LineWidth', 2); 

set(gca,'XTick', 1:iVect, 'XTickLabel', labels)
ylim([min(matrix(:))-.15, max(matrix(:))+.15])



end