mkdir('results');
minE = log2(0.1);
u = {'Esc' 'Dia'};

% ---------------------------------------------------------------------
% load data
% ---------------------------------------------------------------------
load_norm_data = 1;

if (load_norm_data)

    X = importdata('data/Romney_updated_translated_TMM_CPM.txt');
    M = 2.^(X.data);
    Xid = X.textdata(2:end,1);
    cid = X.textdata(1,2:end);
    [~,x] = strtok(cid,'_');
    tm = str2num(cell2mat((strtok(x,'_'))'))';
    S = repmat({'romney'},size(tm));
    [Xid,~,logE] = normalize_fpkm(Xid,Xid,M,zeros(size(M)),S,cid,'results',1,1,2.^minE);
    save('results/data.norm.mat','Xid','logE','cid','tm');

    close all;
end

% ---------------------------------------------------------------------
% fit models
% ---------------------------------------------------------------------
fit_models = 1;

if (fit_models)
    load('results/data.norm.mat','Xid','logE','tm','cid');
    S = strtok(cid,'_');

    % fit a model to each condition (two models)
    for i = 1:max(size(u))
        j = strcmp(S,u{i});
        logEi = logE(:,j);
        tmi = tm(:,j);
        fit_temporal_models(Xid,logEi,tmi,['results.' u{i}],2);
    end

    % fit a model to both conditions
    fit_temporal_models(Xid,logE,tm,'results.two',2);
end

% ---------------------------------------------------------------------
% compare timecourses
% ---------------------------------------------------------------------
compare_tcourses = 1;

if (compare_tcourses)
    load('results/data.norm.mat','Xid','logE','tm','cid');
    compare_timecourses(Xid,logE,tm,u,'two','results',2);
end

