% clear space
clc
clear all
close all

%read data, readtable because file is a table 
data=readtable('Sig_cor_pairs_SNS.csv');

%%%Calculating the line of best fit from survivor data %%%
% number of survivors
n_S = numel(data.x30_days_survival) - nnz(data.x30_days_survival);
n_NS = nnz(data.x30_days_survival);
n = n_S + n_NS;

%Extracting survivor data only for calculating the regression line

results_p1 = NaN(n_S, 2);
results_p2 = NaN(n_S, 2);
results_p3 = NaN(n_S, 2);
results_p4 = NaN(n_S, 2);

index_S = 0;

    for i = 1:162
        if data.x30_days_survival(i) == 0

            index_S = index_S + 1;

            results_p1(index_S,1) = data.wbc(i);
            results_p1(index_S,2) = data.platelets(i);

            results_p2(index_S,1) = data.urea(i);
            results_p2(index_S,2) = data.creatinine(i);

            results_p3(index_S,1) = data.inr(i);
            results_p3(index_S,2) = data.ALT(i);

            results_p4(index_S,1) = data.blood_pH(i);
            results_p4(index_S,2) = data.HCO3(i);

        end
    end


Table_1 = array2table(results_p1,'VariableNames',{'wbc','platelets'});
Table_2 = array2table(results_p2,'VariableNames',{'urea','creatinine'});
Table_3 = array2table(results_p3,'VariableNames',{'inr','ALT'});
Table_4 = array2table(results_p4,'VariableNames',{'blood_pH','HCO3'});

%remove columns with missing rows
Table_pair_1 = rmmissing(Table_1,'DataVariables',{'wbc','platelets'});
Table_pair_2 = rmmissing(Table_2,'DataVariables',{'urea','creatinine'});
Table_pair_3 = rmmissing(Table_3,'DataVariables',{'inr','ALT'});
Table_pair_4 = rmmissing(Table_4,'DataVariables',{'blood_pH','HCO3'});

% calculate m and b for regression line for sig pairs using survivor data
Line_Best_Fit_p1 = linortfit2(Table_pair_1.wbc(:),Table_pair_1.platelets(:));
Line_Best_Fit_p2 = linortfit2(Table_pair_2.urea(:),Table_pair_2.creatinine(:));
Line_Best_Fit_p3 = linortfit2(Table_pair_3.inr(:),Table_pair_3.ALT(:));
Line_Best_Fit_p4 = linortfit2(Table_pair_4.blood_pH(:),Table_pair_4.HCO3(:));

%%% Parenclitic Deviation Calculation %%%
% create table for each pair including both survivors and non survivors
Table_p1_tot = data(:,{'wbc','platelets'});
Table_p2_tot = data(:,{'urea','creatinine'});
Table_p3_tot = data(:,{'inr','ALT'});
Table_p4_tot = data(:,{'blood_pH','HCO3'});

table_list = {Table_p1_tot, Table_p2_tot, Table_p3_tot, Table_p4_tot};

Line_Best_Fit_list = {Line_Best_Fit_p1, Line_Best_Fit_p2, Line_Best_Fit_p3, Line_Best_Fit_p4};

parenclitic_deviations = NaN(n,4);

    for i = 1:4 % numbers of significantly correlated pairs
        table = table2array(table_list{i});
        Line_Best_Fit = Line_Best_Fit_list{i};

        m = Line_Best_Fit(1,1);
        b = Line_Best_Fit(1,2);

        for j = 1:n % numbers of survivors and non survivors

            x = table(j,1);
            y = table(j,2);

            [d] = emilyparencliticdeviation(m,b,x,y);
            parenclitic_deviations(j,i) = d;

        end
    end

%create table that stores subject id, survival status, and parenclitic
%deviations of significantly correlated pairs

parenclitic_deviations_table = array2table(parenclitic_deviations, "VariableNames",{'PD_wbc_platelets','PD_urea_creatinine','PD_inr_ALT','PD_blood_pH_HCO3'});

Final_Table = [data parenclitic_deviations_table];

writetable(Final_Table,'fixed_All_Parenclitic_Dev_S_pairs.csv')




