% Code for finding significantly correlated variables in survivors and 
% nonsurvivors

% Import table data (csv format) - test data used
% Table columns should have contain variable names
clc
clear all
close all

Data_SS = readtable('Sepsis_Survivor_Data.csv'); %Input your own data
Data_NSS = readtable('Sepsis_NonSurvivor_Data.csv'); %Input your own data

z_1 = size(Data_SS,1);
z_2 = size(Data_SS,2);
Data_S = Data_SS(1:z_1,2:z_2);

z_3 = size(Data_NSS,1);
z_4 = size(Data_NSS,2);
Data_NS = Data_NSS(1:z_3,2:z_4);

% Calculate Bonferonni Adjusted P-value
n = size(Data_S,2);
Bonferroni_adjusted_P = 0.05*2/(n*(n-1));

% Calculate Pearsons Correlation Coefficient
[r_S,p_S] = corrcoef(Data_S.Variables,'Rows','pairwise');
[r_NS,p_NS] = corrcoef(Data_NS.Variables,'Rows','pairwise');

% Extract variable names
varNames = Data_S.Properties.VariableNames;

% Make an adjacency matrix with variable names on row and column
P_S_values_Table = array2table(p_S,"RowNames",varNames,"VariableNames",varNames);
P_NS_values_Table = array2table(p_NS,"RowNames",varNames,"VariableNames",varNames);

R_S_values_Table = array2table(r_S,"RowNames",varNames,"VariableNames",varNames);
R_NS_values_Table = array2table(r_NS,"RowNames",varNames,"VariableNames",varNames);

% Create adjacency matrix for survivors and non-survivors, p<Bonferroni_adjusted_P
Adj_matrix_S = R_S_values_Table;
Adj_matrix_NS = R_NS_values_Table;

for r=1:n
    for c=1:n

        if P_S_values_Table{r,c} <= Bonferroni_adjusted_P
             Adj_matrix_S{r,c} = R_S_values_Table{r,c};

        else
            Adj_matrix_S{r,c} = 0;

        end

        if P_NS_values_Table{r,c} <= Bonferroni_adjusted_P
             Adj_matrix_NS{r,c} = R_NS_values_Table{r,c};

        else
            Adj_matrix_NS{r,c} = 0;

        end

    end
end

% Count numbers of correlated variable pairs 
N = nnz(Adj_matrix_S{:,:});
N2 = nnz(Adj_matrix_NS{:,:});

Sig_Pairs_S = NaN(N,2);
Sig_Pairs_NS = NaN(N2,2);

i = 1;
g = 1;

for r = 1:n
   for c = 1:n

        if P_S_values_Table{r,c} <= Bonferroni_adjusted_P
            Sig_Pairs_S(i,1) = r;
            Sig_Pairs_S(i,2) = c;

            i = i + 1;

        else
           
        end

        if P_NS_values_Table{r,c} <= Bonferroni_adjusted_P
            Sig_Pairs_NS(g,1) = r;
            Sig_Pairs_NS(g,2) = c;
            g = g + 1;

        else

        end

   end
end

% Remove Duplicate Pairs in S

Sig_Pairs_S2 = Sig_Pairs_S;

for v=1:N

    if Sig_Pairs_S(v,1) > Sig_Pairs_S(v,2)
        Sig_Pairs_S2(v,2) = Sig_Pairs_S(v,1);
        Sig_Pairs_S2(v,1) = Sig_Pairs_S(v,2);
    else

    end
end

Sig_Pairs_S2_sorted = unique(Sig_Pairs_S2,'rows');

l_1 = size(Sig_Pairs_S2_sorted,1);
l_2 = size(Sig_Pairs_S2_sorted,2);

Sig_Pairs_Survivors = cell(l_1,l_2);

for r = 1:l_1
   for c = 1:l_2
       Sig_Pairs_Survivors{r,c} = varNames{Sig_Pairs_S2_sorted(r,c)};
   end
end

% Remove Duplicate Pairs in NS
Sig_Pairs_NS2 = Sig_Pairs_NS;

for v=1:N2

    if Sig_Pairs_NS(v,1) > Sig_Pairs_NS2(v,2)
        Sig_Pairs_NS2(v,2) = Sig_Pairs_NS(v,1);
        Sig_Pairs_NS2(v,1) = Sig_Pairs_NS(v,2);
    else

    end
end

Sig_Pairs_NS2_sorted = unique(Sig_Pairs_NS2, 'rows');

l_3 = size(Sig_Pairs_NS2_sorted,1);
l_4 = size(Sig_Pairs_NS2_sorted,2);

Sig_Pairs_NonSurvivors = cell(l_3,l_4);

for r = 1:l_3
   for c = 1:l_4
       Sig_Pairs_NonSurvivors{r,c} = varNames{Sig_Pairs_NS2_sorted(r,c)};
   end
end

%Save outputs for significantly correlated pairs in survivors and nonsurvivors
writecell(Sig_Pairs_Survivors,'Significant_Pairs_Survivors.csv');
writecell(Sig_Pairs_NonSurvivors,'Significant_Pairs_Non_Survivors.csv');


