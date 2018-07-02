defaults


nruns =10
all = {};

all_ff = {};
all_ff_err = {};

all_br = {};

all_brturn = {};
all_brnoturn = {};

all_sparse = {};
all_sparse_err = {};

all_pops = {};
all_pops_err = {};


all_ovl = {};
all_ovl_err = {};

all_clusters = {};
all_clusters_err = {};
close all

Conditions = {'10CAP'}
LorG='L';

Turnovers = {0,5,10,15,20};
Patterns = {20};

ninputs=240
%ninputs = 100
npyrs=240


for NCOND=1:length(Conditions)
    CONDITION=Conditions{NCOND}
    
    ss1 = zeros(5, 190);
    
    for NBT=1:length(Turnovers)
        BT=Turnovers{NBT}
        
        
        for NPAT=1:length(Patterns)
            PAT=Patterns{NPAT}
            
            npatterns=PAT
            spks = zeros(npatterns, nruns);
            pops = zeros(npatterns, nruns);
            sparse = zeros(npatterns, nruns);
            
            branch_syns = zeros(npyrs*nbranches, nruns);
            br_hists = zeros(ninputs, 12, nruns);
            clustering = zeros(ninputs, nruns);
            brstrengths = zeros(ninputs, npyrs*nbranches);
            
            brweights = zeros(ninputs, npyrs*nbranches, nruns);
            nrnweights = zeros(ninputs, npyrs, nruns);
            brweightcors = zeros(ninputs, ninputs, nruns);
            brsyncors= zeros(ninputs, ninputs, nruns);
            nrnweightcors = zeros(ninputs, ninputs, nruns);
            
            brcommon = zeros(ninputs, ninputs, nruns);
            
            
            clust_all = {};
            clust_all = cell(9,1);
            
            ISI=120;
            
            bsyns_turn = zeros(npyrs*nbranches, nruns);
            bsyns_noturn = zeros(npyrs*nbranches, nruns);
            
            av = [];

            clusters = [];
            for run = 1:nruns
                
                fn = sprintf('./data/%s%d_0%s_%d_%d/spikesperpattern.dat', CONDITION, PAT, LorG, BT,run-1)
                spk = load( fn);
                
                recallspikes = spk(:, 1:npyrs)/(stimduration/1000);
                pop = recallspikes>CUTOFF; %Hz
                
                spks(:, run) = sum(recallspikes,2);
                
                pops(:,  run) = sum(pop,2);
                

                for npat=1:npatterns-1
                    av(end+1) = sum(pop(npat,:)&pop(npat+1,:));
                end
                
                
                
                for npat=1:npatterns
                    sparse(npat,  run) = trevrolls(recallspikes(npat,:) );

                end
                
                ff = sprintf('./data/%s%d_0%s_%d_%d/synstate.dat', CONDITION, PAT,LorG, BT,run-1);
                ss = load(ff);
                
                for i=1:size(ss,1)
                    bid=ss(i,2);
                    nid=ss(i,3);
                    srcid=ss(i,5);
                    bstrength = ss(i,6);
                    w=ss(i,7);
                    if (srcid >= 0 && bid <= npyrs*nbranches)
                        brweights(srcid+1, bid+1, run) = brweights(srcid+1, bid+1, run) + w;
                        brstrengths(srcid+1, bid+1)=bstrength;
                        nrnweights(srcid+1, nid+1,run) = nrnweights(srcid+1, nid+1,run) + w;
                    end
                    if (srcid >= 0 && bid <= npyrs*nbranches &&  w > 0.7)
                        branch_syns( bid+1, run) = branch_syns( bid+1, run)+1;
                        
                        branchid = mod(bid, nbranches);
                        
                        if (branchid < BT)
                            bsyns_turn(bid+1, run) = bsyns_turn(bid+1, run)+1;
                        else
                            bsyns_noturn(bid+1, run) = bsyns_noturn(bid+1, run)+1;
                        end
                        
                    end
                end
                clusters(end+1) = sum(branch_syns( :, run)>2)/length(branch_syns);
            end
            
            
            all_brturn{ NBT} = bsyns_turn;
            all_brnoturn{ NBT} = bsyns_noturn;
            all_br{NBT} = branch_syns;

            m_p = mean(pops, 2);
            s_p = stderr(pops' )';

            all_pops{NBT} = m_p;
            all_pops_err{NBT} = s_p;

            m_p = mean(spks, 2);
            s_p = stderr(spks' )';

            all_ff{NBT} = m_p/npyrs;
            all_ff_err{NBT} = s_p/npyrs;


            m_p = mean(sparse, 2);
            s_p = stderr(sparse')';

            all_sparse{ NBT} = m_p;
            all_sparse_err{ NBT} = s_p;
            
            
            m_p = mean(av);
            s_p = stderr(av);
            
            ss1(NBT, :) = av';
            
            all_ovl{ NBT} = m_p;
            all_ovl_err{ NBT} = s_p;
            
            m_p = mean(clusters);
            s_p = stderr(clusters);
            
            all_clusters{ NBT} = m_p;
            all_clusters_err{ NBT} = s_p;
        end
        
    end
    
end

% 
% figure()
% z = [];
% 
% for y=1:length(Turnovers)
%     
%     z(:,y) = all_ff{1, y};
% end
% mesh(z)
% tit='FF'
% title(tit);
% export_fig(sprintf('./figs/%s_mesh_%s.pdf',CONDITION, tit), '-transparent')
close all
figure;
hold on
mat  = cell2mat(all_sparse);
mat_err  = cell2mat(all_sparse_err);
errorbar( flipud(mat(:,1)), flipud(mat_err(:,1)), 'r');
errorbar( flipud(mat(:,3)), flipud(mat_err(:,3)), 'b');
errorbar( flipud(mat(:,5)), flipud(mat_err(:,5)), 'g');
legend( 'No turnover', 'Half dendrites', 'All dendrites');
%set(gca, 'XTick', [1:5])
%set(gca, 'XTickLabel', Patterns)
%ylim([.4, .8]);
%title(sprintf('%s', CONDITION));
xlabel('Age of memory (days)');
ylabel('Active Population Sparsity');
%export_fig(sprintf('./figs/%s_%d_SPARSE%s.pdf',CONDITION, PAT, LorG), '-transparent');


figure;
hold on
mat  = cell2mat(all_pops)*100/npyrs;
mat_err  = cell2mat(all_pops_err)*100/npyrs;
errorbar( flipud(mat(:,1)), flipud(mat_err(:,1)), 'r');
errorbar( flipud(mat(:,3)), flipud(mat_err(:,3)), 'b');
errorbar( flipud(mat(:,5)), flipud(mat_err(:,5)), 'g');
legend( 'No turnover', 'Half dendrites', 'All dendrites');
%set(gca, 'XTick', [1:5])
%set(gca, 'XTickLabel', Patterns)
ylim([15, 65]);
%title(sprintf('%s', CONDITION));
xlabel('Age of memory (days)');
ylabel('Active Population (%)');
%export_fig(sprintf('./figs/%s_%d_POPS%s.pdf',CONDITION, PAT,LorG), '-transparent');


figure;
hold on
mat  = cell2mat(all_ff);
mat_err  = cell2mat(all_ff_err);
errorbar( flipud(mat(:,1)), flipud(mat_err(:,1)), 'r');
errorbar( flipud(mat(:,3)), flipud(mat_err(:,3)), 'b');
errorbar( flipud(mat(:,5)), flipud(mat_err(:,5)), 'g');
legend( 'No turnover', 'Half dendrites', 'All dendrites');
%set(gca, 'XTick', [1:5])
%set(gca, 'XTickLabel', Patterns)
%ylim([.4, 1]);
%ylim([5, 40]);
%title(sprintf('%s', CONDITION));
xlabel('Age of memory (days)');
ylabel('Avg Firing Rate (Hz)');

%export_fig(sprintf('./figs/%s_%d_FF%s.pdf',CONDITION, PAT,LorG), '-transparent');



figure;
hold on
mat  = cell2mat(all_clusters);
mat_err  = cell2mat(all_clusters_err);
errorbar( fliplr(mat)*100, fliplr(mat_err)*100, 'b');
hold off;
	set(gca, 'XTick', [1:5])
	set(gca, 'XTickLabel', [0:5:20])
%ylim([.4, 1]);
%%title(sprintf('%s', CONDITION));
xlabel('Number of Dendrites with turnover');
ylabel('Branches with clusters (%)');
%export_fig(sprintf('./figs/%s_%d_CLUST%s.pdf',CONDITION, PAT,LorG), '-transparent');
csvwrite(sprintf('./figs/%s_%d_CLUST%s.csv',CONDITION, PAT,LorG), [fliplr(mat)*100; fliplr(mat_err)*100] )


figure;
hold on
mat  = cell2mat(all_ovl)*100/npyrs;
mat_err  = cell2mat(all_ovl_err)*100/npyrs;
errorbar( (mat), (mat_err), 'b');
hold off;
	set(gca, 'XTick', [1:5])
	set(gca, 'XTickLabel', [0:5:20])
%ylim([.4, 1]);
%%title(sprintf('%s', CONDITION));
xlabel('Number of Dendrites with turnover');
ylabel('Avg population overlap (%)');
%export_fig(sprintf('./figs/%s_%d_OVL%s.pdf',CONDITION, PAT,LorG), '-transparent');
csvwrite(sprintf('./figs/%s_%d_OVL%s.csv',CONDITION, PAT,LorG), [mat; mat_err] )




mat  = cell2mat(all_sparse);
csvwrite(sprintf('./figs/%s_%d_SPAR%s.csv',CONDITION, PAT,LorG), mat)

%errorbar( mean(mat), stderr(mat), 'r');

%legend( 'No turnover', 'Half dendrites', 'All dendrites');
%set(gca, 'XTick', [1:5])
%set(gca, 'XTickLabel', Patterns)
%ylim([.4, .8]);
%title(sprintf('%s', CONDITION));
%xlabel('Age of memory (days)');
%ylabel('Active Population Sparsity');
%export_fig(sprintf('./figs/%s_SPARSE.pdf',CONDITION), '-transparent');


