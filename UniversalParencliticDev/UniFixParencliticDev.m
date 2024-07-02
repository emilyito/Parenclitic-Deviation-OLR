% Code for calculating Parenclitic Network for any data in a table
% Requires a cell array containing information about significantly
% correlated variable pairs (Use -> FindSigCorPairsSNS.m)

clc
clear all
close all

Data_S = readtable('Sepsis_Survivor_Data.csv'); %Input your own data
Data_NS = readtable('Sepsis_NonSurvivor_Data.csv'); %Input your own data


l = size(Data_S,1);
l2 = size(Data_NS,1);

%%% Load the datas for significant pairs in survivors and nonsurvivors %%%
Data_Survivors_Sig_Cell = readcell("Significant_Pairs_Survivors.csv"); 
Data_NonSurvivors_Sig_Cell = readcell("Significant_Pairs_Non_Survivors.csv");

%%% Calculate PD for significant survivor/nonsurvivor pairs 
Sig_Pairs_List = {Data_Survivors_Sig_Cell, Data_NonSurvivors_Sig_Cell};
Sig_Pair_Group = {'Sig_S_Pair', 'Sig_NS_Pair'};

for index = 1:2
    
    z = size(Sig_Pairs_List{index},1);

    P_Dev_All_S = NaN(l,z); %Significant survivor pair compare with all survivors
    P_Dev_All_NS = NaN(l2,z); %Significant survivor pair compare with all nonsurvivors
    PD_Variable_Names = cell(z,1);

    Data_Sig_Cell = Sig_Pairs_List{index};

    for  p = 1:z

%%% Extract data for significantly correlated variables from the original table
        x = Data_Sig_Cell(p,1); % Extract data on significantly correlated variables
        y = Data_Sig_Cell(p,2);

        a = string(x); % convert to string to extract data from the original table
        b = string(y);

        Data_XY_S = table; % Open table for all survivors
        Data_XY_NS = table; % Open table for all nonsurvivors

        Data_XY_S{1:l,a} = Data_S.(a); % survivors - Variable x
        Data_XY_S{1:l,b} = Data_S.(b); % survivors - Variable y

        Data_XY_NS{1:l2,a} = Data_NS.(a); % nonsurvivors - Variable x
        Data_XY_NS{1:l2,b} = Data_NS.(b); % nonsurvivors - Variable y

%%% Clean data by removing rows with NaN for survivor data (reference population)
        subject_ID_S = Data_S{:,1};
        subject_ID_S_Table = array2table(subject_ID_S);

        subject_ID_NS = Data_NS{:,1};
        subject_ID_NS_Table = array2table(subject_ID_NS);

        Data_XY_S = [subject_ID_S_Table Data_XY_S];
        Data_XY_NS = [subject_ID_NS_Table Data_XY_NS];
    
        Data_XY_S_clean = rmmissing(Data_XY_S);

        %calculate the size of the cleaned survivor data with NaN rows removed
        n=size(Data_XY_S_clean,1);

        %%% Calculating the line of best fit (LOBF) from survivor data %%%
        % Use orthogonal linear regression
        Line_Best_Fit = linortfit2(Data_XY_S_clean.(a),Data_XY_S_clean.(b));

        %%% Parenclitic Deviation Calculation %%%
        m = Line_Best_Fit(1,1);
        c = Line_Best_Fit(1,2);

        parenclitic_deviations_S = NaN(l,1);
        parenclitic_deviations_NS = NaN(l2,1);
    
        % Calculate PD for Survivor Data
        for i = 1:l

            x_S = Data_XY_S.(a)(i);
            y_S = Data_XY_S.(b)(i);

           [d] = emilyparencliticdeviation(m,c,x_S,y_S);
           parenclitic_deviations_S(i,1) = d;

        end

    % Calculate PD for Nonsurvivor Data
        for i = 1:l2

            x_NS = Data_XY_NS.(a)(i);
            y_NS = Data_XY_NS.(b)(i);

           [d] = emilyparencliticdeviation(m,c,x_NS,y_NS);
           parenclitic_deviations_NS(i,1) = d;

        end

        P_Dev_All_S(1:l,p) = parenclitic_deviations_S(1:l,1);
        P_Dev_All_NS(1:l2,p) = parenclitic_deviations_NS(1:l2,1);

        % Add New Variable Names on to the list
        PD_ab = strcat('PD_', a, '_', b);
        PD_Variable_Names(p,1) = cellstr(PD_ab);

    end

    %%% Convert PD data in matrix to table
    % variable name automatically generated using variable pair data

    P_Dev_All_S_Table = array2table(P_Dev_All_S);
    P_Dev_All_NS_Table = array2table(P_Dev_All_NS);

    for  p = 1:z

       P_Dev_All_S_Table.Properties.VariableNames{p} = char(PD_Variable_Names(p,1));
       P_Dev_All_NS_Table.Properties.VariableNames{p} = char(PD_Variable_Names(p,1));

    end

    P_Dev_All_S_Table = [subject_ID_S_Table P_Dev_All_S_Table];
    P_Dev_All_NS_Table = [subject_ID_NS_Table P_Dev_All_NS_Table];

    newTableName_S = string(strcat('P_Dev_All_', Sig_Pair_Group{index}, '_S_Table'));
    writetable(P_Dev_All_S_Table, string(strcat(newTableName_S,'.csv')));

    newTableName_NS = string(strcat('P_Dev_All_', Sig_Pair_Group{index}, '_NS_Table'));
    writetable(P_Dev_All_NS_Table, string(strcat(newTableName_NS,'.csv')));

end

