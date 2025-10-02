%folder = '/Volumes/Caroline_backup/Ox68/';
folder = strcat(plot_dir,'/tables/');
all_stems = ["hr_mean","hr_std","sats_mean","sats_std","rr_mean","rr_std"];
for fname_stem = all_stems
    %fname_stem = "sats_mean";
    
    %data_table=readtable(strcat(folder,"hr_mean.xlsx"));
    %data_table=readtable(strcat(folder,"hr_std.xlsx"));
 data_table = readtable(strcat(folder,fname_stem,sprintf('_%s.xlsx',dataset_label)));
    %data_table=readtable(strcat(folder,"hr_mean.xlsx"));
    %data_table=readtable(strcat(folder,"hr_std.xlsx"));
    %data_table=readtable(strcat(folder,"sats_mean.xlsx"));
    %data_table=readtable(strcat(folder,"rr_mean.xlsx"));
    %data_table=readtable(strcat(folder,"sats_std.xlsx"));
    %data_table=readtable(strcat(folder,"rr_std.xlsx"));
    
    %demographics
    sex=categorical(data_table.BIRTH_Gender);
    PMA=data_table.PMA;
    weight=data_table.TESTOD_MostRecentWeight;
    ventilation=categorical(data_table.TESTOD_CurrentVentilation);
    
    %transfusion
    starthb=data_table.StartHB;
    endhb=data_table.EndHB;
    diffhb=endhb-starthb;%diffhb=data_table.diffHb;
    PNA=data_table.PNA;
    
    %metrics
    increase=data_table.AverageSignalIncrease;
    decrease=data_table.AverageSignalDecrease;
    overall=increase-decrease;
    
    %%
    % t=table(PMA,weight,ventilation,diffhb,increase);
    % mdl=fitlm(t,'ResponseVar','increase')
    % 
    % t=table(PMA,weight,ventilation,diffhb,decrease);
    % mdl=fitlm(t,'ResponseVar','decrease')
    
    response=overall;
    % response=-decrease;
    % response=increase;
    
    %t=table(weight,ventilation,starthb,response);
    %t=table(weight,ventilation,diffhb,starthb,response);
    t=table(diffhb,PNA,PMA, response);
    %mdl=fitlm(t,'response~ventilation*weight*diffhb')
    %mdl=fitlm(t,'response~ventilation+weight+diffhb');
    mdl=fitlm(t,'response~PNA+PMA+diffhb');
    

    figure; plotAdjustedResponse(mdl,'PNA')
    
    %%
    
    figure; 
    subplot(3,3,1); scatter(sex,response,'filled'); xlabel('sex','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,2); scatter(PMA,response,'filled'); lsline; xlabel('PMA','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,3); scatter(ventilation,response,'filled'); xlabel('ventilation mode','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,4); scatter(weight,response,'filled'); lsline; xlabel('weight','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,5); scatter(PNA,response,'filled'); lsline; xlabel('PNA','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,6); scatter(starthb,response,'filled'); lsline; xlabel('starthb','fontsize',15); set(gca,'fontsize',15)
    subplot(3,3,7); scatter(diffhb,response,'filled'); lsline; xlabel('diffhb','fontsize',15); set(gca,'fontsize',15)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    print(strcat(folder,fname_stem,'_correlations.pdf'),'-dpdf','-bestfit')
    close
end