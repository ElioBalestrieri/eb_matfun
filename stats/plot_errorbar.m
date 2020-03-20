function plot_errorbar(matrix, labels, rndmlevel)

figure; hold on;

avg_ = mean(matrix);
% bar(avg_, 'FaceColor', [100 192 238]/255); hold on
err_ = std(matrix)/sqrt(size(matrix,1));

for iVect = 1:size(matrix, 2)
    this_vect = matrix(:,iVect);
    xsc_ = (randn(numel(this_vect),1))/30 +iVect;
    scatter(xsc_, this_vect, 20, [100 192 238]/255, 'filled')
    
    [~, p, ~, stat] = ttest(this_vect, rndmlevel, 'Tail', 'right');
    
    stringsign1 = sprintf('t = %4.4f', stat.tstat);
    stringsign2 = sprintf('p = %4.4f', p);

    text(iVect-.15, max(matrix(:))+.1, {stringsign1, stringsign2})
    
end

errorbar(1:length(avg_), avg_, err_, '.k', 'LineWidth', 2); 
plot(0:iVect+1, ones(iVect+2,1)*rndmlevel,'r', 'LineWidth', 2)


set(gca,'XTick', 1:iVect, 'XTickLabel', labels)
ylim([min(matrix(:))-.15, max(matrix(:))+.15])




end