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

for maxShell=1:14

    M = maxShell*(maxShell+1)/2;
    diffEq = strcat('HKC',num2str(M));
    name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_SmallRay.mat');
    allData = load(name);
    surf(allData.RoGrid,allData.RaGrid,allData.heatTransportGrid,'FaceColor',colorScheme(mod(maxShell-1,7)+1,:),'FaceAlpha',0.75); hold on;
    legendLabels{maxShell} = num2str(M);
end

legend(legendLabels,'NumColumns',2);

for i =1:7
    fig = figure; hold on;
    for maxShell=1:14
        M = maxShell*(maxShell+1)/2;
        diffEq = strcat('HKC',num2str(M));
        name = strcat('../SavedData/HeatTransport/',diffEq,'_HeatTransport_SmallRay.mat');
        allData = load(name);
        if(maxShell<14)
            plot(allData.RaGrid(:,i),allData.heatTransportGrid(:,i),'LineWidth',1.5); %'Color',colorScheme(2,:)
        else
            plot(allData.RaGrid(:,i),allData.heatTransportGrid(:,i),'Color',colorScheme(4,:),'LineWidth',3.5);
            plot(allData.RaGrid(:,i),0.9*allData.heatTransportGrid(:,i),'Color',colorScheme(6,:),'LineWidth',3.5);
            plot(allData.RaGrid(:,i),1.1*allData.heatTransportGrid(:,i),'Color',colorScheme(6,:),'LineWidth',3.5);
        end
    end
    legend(legendLabels,'NumColumns',3);
    xlim([0 5000]); ylim([.85 5.25]);
    title(strcat('Rot = ',num2str((i-1)*50))); fontsize(fig, 16, "points");
end
