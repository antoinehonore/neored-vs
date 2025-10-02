%folder = '/Volumes/Caroline_backup/Ox68/';
folder = strcat(plot_dir,'/tables/');
all_stems = ["hr_mean","hr_std","sats_mean","sats_std","rr_mean","rr_std"];
for fname_stem = all_stems
    %fname_stem = "sats_mean";
    
    %data_table=readtable(strcat(folder,"hr_mean.xlsx"));
    %data_table=readtable(strcat(folder,"hr_std.xlsx"));
    data_table = readtable(strcat(folder,fname_stem,sprintf('_%s.xlsx',dataset_label)));
    %data_table=readtable(strcat(folder,"rr_mean.xlsx"));
    %data_table=readtable(strcat(folder,"sats_std.xlsx"));
    %data_table=readtable(strcat(folder,"rr_std.xlsx"));
    %data_table.BIRTH_Gender = changem(data_table.BIRTH_Gender,77,'M');
    %demographics
    sex=categorical(data_table.BIRTH_Gender);
    PMA=data_table.PMA/7;
    PNA=data_table.PNA;
    weight=data_table.TESTOD_MostRecentWeight;
    ventilation=categorical(data_table.TESTOD_CurrentVentilation);
    
    %transfusion
    starthb=data_table.StartHB;
    endhb=data_table.EndHB;
    diffhb=endhb-starthb;
    
    
    %metrics
    increase_responder=data_table.IncreaseResponse;
    decrease_responder=data_table.DecreaseResponse;
    responder=increase_responder+decrease_responder;
    
    response=zeros(length(increase_responder),1);
    
    response(find(responder==0))=1;
    response(find(decrease_responder))=2;
    response(find(increase_responder))=3;
    
    %%

    fig=figure; boxplot(PMA,response); ylabel('PMA (weeks)','fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'}); set(gca,'fontsize',15)
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_pma.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_pma.jpg'),'-dpdf','-bestfit')

    fig=figure; boxplot(PNA,response); ylabel('PNA (days)','fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'}); set(gca,'fontsize',15)
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_pna.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_pma.jpg'),'-dpdf','-bestfit')


    fig=figure; boxplot(weight,response); ylabel('Weight (kg)','fontsize',15); set(gca,'fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'});
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    orient(fig, 'landscape')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_weight.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_weight.jpg'),'-dpdf','-bestfit')

    fig=figure; boxplot(starthb,response); ylabel('Start Hb (g/l)','fontsize',15); set(gca,'fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'});
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    orient(fig, 'landscape')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_starthb.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_starthb.jpg'), '-dpdf','-bestfit')

    fig=figure; boxplot(endhb,response); ylabel('End Hb (g/l)','fontsize',15); set(gca,'fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'});
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    orient(fig, 'landscape')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_endhb.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_endhb.jpg'),'-dpdf','-bestfit')

    fig=figure; boxplot(diffhb,response); ylabel('Increment Hb (g/l)','fontsize',15); set(gca,'fontsize',15); set(gca,'XTickLabel',{'Non','Decrease','Increase'});
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    title(fname_stem, 'Interpreter', 'none')
    orient(fig, 'landscape')
    exportgraphics(gcf, strcat(folder,fname_stem,'_boxplots_incrhb.jpg'));
    %print(strcat(folder,fname_stem,'_boxplots_incrhb.jpg'),'-dpdf','-bestfit')
    close all
end
%%
%!pdftk ../results/tables/*pma.jpg cat output ../results/tables/pma_all.jpg
%!pdftk ../results/tables/*weight.jpg cat output ../results/tables/weight_all.jpg
%!pdftk ../results/tables/*starthb.jpg cat output ../results/tables/starthb_all.jpg
%!pdftk ../results/tables/*endhb.jpg cat output ../results/tables/endhb_all.jpg
%!pdftk ../results/tables/*incrhb.jpg cat output ../results/tables/incrhb_all.jpg
