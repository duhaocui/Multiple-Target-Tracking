%% Load simulation data

filePattern = fullfile(pwd,'*.mat');
simulationData = dir(filePattern);
for k = 1:length(simulationData)
    baseFileName = simulationData(k).name;
    baseFileName = baseFileName(1:end-4);
    data = load(baseFileName);
    v = genvarname(baseFileName, who);
    eval([v '= data.averGospa;']);
end
%%
x = '_30_75';
glmb = eval(strcat('glmb',x));
lmb = eval(strcat('lmb_merge',x));
pmbm = eval(strcat('pmbm',x));
pmbm_recycle = eval(strcat('pmbm_recycle',x));
pmb = eval(strcat('pmb',x));
pmb_recycle = eval(strcat('pmb_recycle',x));

figure(1)
for i = 1:4
    subplot(2,2,i)
    plot(glmb(:,i),'LineWidth',1);
    hold on
    grid on
    box on
    plot(lmb(:,i),'LineWidth',1);
    plot(pmbm(:,i),'LineWidth',1);
    plot(pmbm_recycle(:,i),'LineWidth',1);
    plot(pmb(:,i),'LineWidth',1);
    plot(pmb_recycle(:,i),'LineWidth',1);
    lgd = legend('\delta-GLMB','LMB','PMBM w/o recycling','PMBM w/ recycling','PMB w/o recycling','PMB w/ recycling','Location','best');
    lgd.FontSize = 6;
    xlabel('time step')
    switch i
        case 1
            ylabel('Mean GOSPA (Total)');
        case 2
            ylabel('Mean GOSPA (Loc)')
        case 3
            ylabel('Mean GOSPA (Missed)')
        case 4
            ylabel('Mean GOSPA (False)')
    end
end

