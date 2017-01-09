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
x = '_10_75';
glmb = eval(strcat('glmb',x));
lmb = eval(strcat('lmb_merge',x));
pmbm = eval(strcat('pmbm',x));
pmbm_recycle = eval(strcat('pmbm_recycle',x));
pmb = eval(strcat('pmb',x));
pmb_recycle = eval(strcat('pmb_recycle',x));

figure(1)
for i = 1:4
    subplot(2,2,i)
    hold on
    grid on
    xlabel('time step')
    ylabel('Mean GOSPA')
    plot(glmb(:,i),'LineWidth',1);
    plot(lmb(:,i),'LineWidth',1);
    plot(pmbm(:,i),'LineWidth',1);
    plot(pmbm_recycle(:,i),'LineWidth',1);
    plot(pmb(:,i),'LineWidth',1);
    plot(pmb_recycle(:,i),'LineWidth',1);
    legend('GLMB','LMB','PMBM w/o recycling','PMBM w/ recycling','PMB w/o recycling','PMB w/ recycling')
end
%%
% N = truth.total_tracks;
% tracks = cell(N,1);
% for i = 1:N
%     tracks{i} = zeros(2,101);
%     for j = 1:101
%         tracks{i}(:,j) = truth.X{j}([1,3],i);
%     end
% end
% figure(1)
% subplot(2,1,1)
% hold on
% for i = 1:N
%     plot(0:100,tracks{i}(1,:));
% end
% xlabel('time step')
% ylabel('Position(x)')
% subplot(2,1,2)
% hold on
% for i = 1:N
%     plot(0:100,tracks{i}(2,:));
% end
% xlabel('time step')
% ylabel('Position(y)')
