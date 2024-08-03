% Code for visualising the parenclitic network

clc
clear all
close all

% Read table to calculate numbers of variables

Data_SS = readtable('Sepsis_Survivor_Data.csv'); % input your own data
z_1 = size(Data_SS,2);
Data_S = Data_SS(:,2:z_1);
z_2 = size(Data_S,2);

% Read data for significant pairs

Data_Survivors_Sig_Cell = readcell("Significant_Pairs_Survivors.csv"); 
Data_NonSurvivors_Sig_Cell = readcell("Significant_Pairs_Non_Survivors.csv");

SNS_Sig_Cell = [Data_Survivors_Sig_Cell; Data_NonSurvivors_Sig_Cell];
SNS_Sig_Cell_table = cell2table(SNS_Sig_Cell);
SNS_Sig_Cell_table = unique(SNS_Sig_Cell_table(:,1:2),'rows');
SNS_Sig_Cell = table2cell(SNS_Sig_Cell_table);

n = size(SNS_Sig_Cell,1);

% Read data for the parenclitic deviation computed from
% UniFixParencliticDev.m

PD_Data_Sig_S_Pair_S = readtable('P_Dev_All_Sig_S_Pair_S_Table.csv'); 
PD_Data_Sig_NS_Pair_S = readtable('P_Dev_All_Sig_NS_Pair_S_Table.csv');
PD_Data_Sig_S_Pair_NS = readtable('P_Dev_All_Sig_S_Pair_NS_Table.csv');
PD_Data_Sig_NS_Pair_NS = readtable('P_Dev_All_Sig_NS_Pair_NS_Table.csv');

% combine data for different patient groups into one cell/table 

n_1 = size(PD_Data_Sig_S_Pair_S,2);
n_2 = size(PD_Data_Sig_NS_Pair_S,2);
n_3 = size(PD_Data_Sig_S_Pair_NS,2);
n_4 = size(PD_Data_Sig_NS_Pair_NS,2);

subject_ID_S = PD_Data_Sig_S_Pair_S{:,1};
subject_ID_S_Table = array2table(subject_ID_S);

subject_ID_NS = PD_Data_Sig_S_Pair_NS{:,1};
subject_ID_NS_Table = array2table(subject_ID_NS);

variableNames_Sig_S_Pair_S = PD_Data_Sig_S_Pair_S.Properties.VariableNames;
variableNames_Sig_NS_Pair_S = PD_Data_Sig_NS_Pair_S.Properties.VariableNames;

common_variables_S = intersect(variableNames_Sig_S_Pair_S, variableNames_Sig_NS_Pair_S);
common_variables_S2 = string(common_variables_S);

PD_Data_Sig_NS_Pair_S = removevars(PD_Data_Sig_NS_Pair_S, common_variables_S2);

variableNames_Sig_S_Pair_NS = PD_Data_Sig_S_Pair_NS.Properties.VariableNames;
variableNames_Sig_NS_Pair_NS = PD_Data_Sig_NS_Pair_NS.Properties.VariableNames;

common_variables_NS = intersect(variableNames_Sig_S_Pair_NS, variableNames_Sig_NS_Pair_NS);
common_variables_NS2 = string(common_variables_NS);

PD_Data_Sig_NS_Pair_NS = removevars(PD_Data_Sig_NS_Pair_NS, common_variables_NS2);

P_Dev_All_S_Table = [PD_Data_Sig_S_Pair_S PD_Data_Sig_NS_Pair_S];
P_Dev_All_NS_Table = [PD_Data_Sig_S_Pair_NS PD_Data_Sig_NS_Pair_NS];

l = size(P_Dev_All_S_Table,2);

% Make an adjacency matrix with variable names on row and column
% no correlation = 0
% yes correlation = 1

% Extract variable names
varNames = Data_S.Properties.VariableNames;

Adjacency_matrix_S = zeros(z_2,z_2);
Adjacency_matrix_NS = zeros(z_2,z_2);

Adj_S_Table = array2table(Adjacency_matrix_S,"RowNames",varNames,"VariableNames",varNames);
Adj_NS_Table = array2table(Adjacency_matrix_NS,"RowNames",varNames,"VariableNames",varNames);

for i = 1:n

    a = string(SNS_Sig_Cell(i,1));
    b = string(SNS_Sig_Cell(i,2));

    Adj_S_Table{a,b} = 1;
    Adj_NS_Table{a,b} = 1;

    Adj_S_Table{b,a} = 1; % this is because adj_mat should be symmetric
    Adj_NS_Table{b,a} = 1;

end

% Find Average PD values for survivors and nonsurvivors

% omitnan calculates the mean after removing NaN
func = @(x) mean(x,"omitnan");

Mean_P_Dev_All_S_Table = varfun(func,P_Dev_All_S_Table(:,2:l));
Mean_P_Dev_All_NS_Table = varfun(func,P_Dev_All_NS_Table(:,2:l));

variableNames_All_Pairs = Mean_P_Dev_All_S_Table.Properties.VariableNames;

Adj_S_Table_wPD = Adj_S_Table;
Adj_NS_Table_wPD = Adj_NS_Table;

for z = 1:n

    str = variableNames_All_Pairs{z};
    expression = 'Fun_PD_(.*)_(.*)';
    tokens = regexp(str, expression, 'tokens');

    a = tokens{1}{1};
    b = tokens{1}{2};

    a = string(a);
    b = string(b);

    Adj_S_Table_wPD{a,b} = Mean_P_Dev_All_S_Table{1,z};
    Adj_NS_Table_wPD{a,b} = Mean_P_Dev_All_NS_Table{1,z};

    Adj_S_Table_wPD{b,a} = Mean_P_Dev_All_S_Table{1,z};
    Adj_NS_Table_wPD{b,a} = Mean_P_Dev_All_NS_Table{1,z};

end

%convert adjacency matrix into graphs

Adj_S_Array_wPD = table2array(Adj_S_Table_wPD);

Adj_NS_Array_wPD = table2array(Adj_NS_Table_wPD);

G_S = graph(Adj_S_Array_wPD);

l=5*abs(G_S.Edges.Weight);

subplot(1,2,1)
Network_S = plot(G_S,'EdgeLabel',G_S.Edges.Weight,'Layout','circle', 'LineWidth',l, EdgeFontSize=10, NodeFontSize=12);

for i = 1:z_2
labelnode(Network_S,i,varNames{i});
end 

title('Survivors',FontSize=15)

axis off;

G_NS=graph(Adj_NS_Array_wPD);

l2=5*abs(G_NS.Edges.Weight);

subplot(1,2,2)
Network_NS = plot(G_NS,'EdgeLabel',G_NS.Edges.Weight,'Layout','circle','LineWidth',l2, EdgeFontSize=10, NodeFontSize=12);

for i = 1:z_2
labelnode(Network_NS,i,varNames{i});
end 

title('Non-Survivors', FontSize=15)

axis off;
