% Plotting heat transport across different HKC models
% by Roland Welter

clear all, close all, clc

colorScheme(1,:)=[0.4940 0.1840 0.5560];
colorScheme(2,:)=[0 0.4470 0.7410];
colorScheme(3,:)=[0.3010 0.7450 0.9330];
colorScheme(4,:)=[0.4660 0.6740 0.1880];
colorScheme(5,:)=[0.9290 0.6940 0.1250];
colorScheme(6,:)=[0.8500 0.3250 0.0980];
colorScheme(7,:)=[0.6350 0.0780 0.1840];


for maxShell=8:14

    M = maxShell*(maxShell+1)/2;
    diffEq = strcat('HKC',num2str(M));
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport.mat');
    allData = load(name);
    surf(allData.RoGrid,allData.RaGrid,allData.heatTransportGrid,'FaceColor',colorScheme(maxShell-7,:),'FaceAlpha',0.75); hold on;
    legendLabels{maxShell-7} = strcat('HKC ',num2str(M));
end

legend(legendLabels); 

for i =1:11
    fig = figure; hold on;
    for maxShell=8:14
        M = maxShell*(maxShell+1)/2;
        diffEq = strcat('HKC',num2str(M));
        name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport.mat');
        allData = load(name);
        plot(allData.RaGrid(:,i),allData.heatTransportGrid(:,i),'Color',colorScheme(maxShell-7,:),'LineWidth',2.5);
    end
    legend(legendLabels);
    xlim([0 25*10^3]);
    title(strcat('Rot = ',num2str((i-1)*10^2))); fontsize(fig, 16, "points");
end
