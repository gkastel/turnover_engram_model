
figure;
hold on

mat = load('./figs/10CAP_20_OVLG.csv')
errorbar( mat(1,:), mat(2,:), 'g');

mat = load('./figs/10CAP_20_OVLL.csv')
errorbar( mat(1,:), mat(2,:), 'r');

legend('Somatic PRP synthesis', 'Dendritic PRP synthesis');
hold off;
	set(gca, 'XTick', [1:5])
	set(gca, 'XTickLabel', [0:5:20])
%ylim([.4, 1]);
%%title(sprintf('%s', CONDITION));
ylim([18,26]);
xlabel('Number of Dendrites with turnover');
ylabel('Average population overlap (%)');
export_fig(sprintf('./figs/%s_OVLBOTH.pdf',CONDITION), '-transparent');


figure;
hold on

mat = load('./figs/10CAP_20_SPARG.csv');
%mat = mat(19:20, :);
mean(mat)
stderr(mat)

errorbar( mean(mat), stderr(mat), 'g');

mat = load('./figs/10CAP_20_SPARL.csv');
mean(mat)
stderr(mat)

%mat = mat(19:20, :);
errorbar( mean(mat), stderr(mat), 'r');

legend('Somatic PRP synthesis', 'Dendritic PRP synthesis');
hold off;
set(gca, 'XTick', [1:5])
set(gca, 'XTickLabel', [0:5:20])
ylim([.55, .8]);
%%title(sprintf('%s', CONDITION));
%ylim([16,28]);
xlabel('Number of Dendrites with turnover');
ylabel('Average population sparsity (%)');
export_fig(sprintf('./figs/%s_SPARSBOTH.pdf',CONDITION), '-transparent');
%csvwrite(sprintf('./figs/%s_%d_OVL%s.csv',CONDITION, PAT,LorG), [mat; mat_err] )

n = zeros(npyrs);
hi  = 0.3 * npyrs
n(1:hi) = 1;
com = [];
for i=1:1000;
    
    ix = randperm(npyrs);
    n1 = n(ix);
    ix = randperm(npyrs);
    n2 = n(ix);
    com(end+1) = sum(n1 & n2);
end;
100*mean(com)/npyrs



