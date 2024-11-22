%% 1. normalization of protein/phosphopeptide data for clustering analysis
% 1) import raw abundance data as 'allo_raw'
allo_raw(allo_raw==0)=NaN;

% 2) normalized between TMT channels
allo_ctl=allo_raw(:,1:18);
allo_tx=allo_raw(:,19:36);

tm_ctl_q=quantilenorm(allo_ctl);
tm_tx_q=quantilenorm(allo_tx);

% 3) log2-transformation
tm_log_ctl=log2(tm_ctl_q+1);
tm_log_tx=log2(tm_tx_q+1);

% 4) fold-change calculation relative to global reference
tm_ratio_ctl=tm_log_ctl-repmat(tm_log_ctl(:,1),1,18);tm_ratio_ctl=tm_ratio_ctl(:,2:end); 
tm_ratio_tx=tm_log_tx-repmat(tm_log_tx(:,1),1,18);tm_ratio_tx=tm_ratio_tx(:,2:end);

% 5) normalize between TMT sets
tm_rs_ctl=reshape(tm_ratio_ctl,[size(tm_ctl_q,1)*17,1]);
tm_rs_tx=reshape(tm_ratio_tx,[size(tm_ctl_q,1)*17,1]);

tm_m_q=quantilenorm([tm_rs_ctl,tm_rs_tx]);

allo_norm=[reshape(tm_m_q(:,1),size(tm_ctl_q,1),17),reshape(tm_m_q(:,2),size(tm_ctl_q,1),17)];
allo_norm=quantilenorm(allo_norm);

%% 2. Treated-Control fold-change calculation and expression filtering
FC=allo_norm(:,18:34)-allo_norm(:,1:17);
FC_q=quantilenorm(FC);

tm_ind_exp=sum(isnan(FC_q),2)==0;
FC_q=FC_q(tm_ind_exp,:);

%% 3. ONMF clustering
% 1) Filter Top MAD protein/phosphopeptides
FC_q_mad=mad(FC_q,1,2);
ind_mad{1,1}=find(FC_q_mad>prctile(FC_q_mad,90));
ind_mad{2,1}=find(FC_q_mad>prctile(FC_q_mad,80));
ind_mad{3,1}=find(FC_q_mad>prctile(FC_q_mad,70));

% 2) ONMF clustering
for i=1:3
    i
    tm_data=FC_q(ind_mad{i,1},:);
    for j=2:8
        j
        [coph_cor{i}(j-1),ave_C{i}{j-1},~,~]=aoNMF_subtyping(auto(tm_data'),j,1000,100);
    end
end

% 3) Cluster consensus value matrix
CG_mad30_c2=clustergram(ave_C{1,3}{1,1},'cluster',3,'standardize',3,'Linkage','complete');



