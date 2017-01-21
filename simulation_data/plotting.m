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
x = '_coal_10_98';
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
%%
N = truth.total_tracks;
tracks = cell(N,1);
for i = 1:N
    tracks{i} = zeros(2,101);
    for j = 1:101
        tracks{i}(:,j) = truth.X{j}([1,3],i);
    end
end
figure(1)
subplot(2,1,1)
hold on
for i = 1:N
    plot(0:100,tracks{i}(1,:));
end
xlabel('time step')
ylabel('Position(x)')
subplot(2,1,2)
hold on
for i = 1:N
    plot(0:100,tracks{i}(2,:));
end
xlabel('time step')
ylabel('Position(y)')

% for i = 1:N
%     plot3(tracks{i}(1,:),tracks{i}(2,:),0:100);
%     hold on
% end

%%
[X_track,k_birth,k_death]= extract_tracks(truth.X,truth.track_list,truth.total_tracks);

%plot ground truths
limit= [ model.range_c(1,1) model.range_c(1,2) model.range_c(2,1) model.range_c(2,2) 1 100];
figure; truths= gcf;
for i=1:truth.total_tracks
    Pt= X_track(:,k_birth(i):1:k_death(i),i); Pt=Pt([1 3],:);
    plot3( Pt(1,:),Pt(2,:),k_birth(i):k_death(i)); 
    hold on
end
axis(limit); title('Ground Truths');
grid on
