clc; clear;

N = input('Enter number of nodes: ');
n = input('Enter number of trials per edge: ');
prior_H1 = input('Enter prior probability of edge being present (0-1) for MAP: ');

edges = nchoosek(1:N,2);
numedges = size(edges,1);
edgenames = cellstr(strcat(string(edges(:,1)),'-',string(edges(:,2))));

fprintf('\nEdges in the network:\n');
disp(edges);

fprintf('Enter observed trials for each edge (0=absent, 1=present) as a row vector of length %d:\n', n);
datamatrix = zeros(numedges, n);
for i = 1:numedges
    datamatrix(i,:) = input(sprintf('Edge %s: ', edgenames{i}));
end

dataTable = array2table(datamatrix,'VariableNames',strcat('Trial', string(1:n)),'RowNames',edgenames);
disp('Observed Data Table:');
disp(dataTable);

X = sum(datamatrix,2);
p1 = input('Enter P(edge works|H1) (high, e.g., 0.8): ');
p0 = input('Enter P(edge works|H0) (low, e.g., 0.2): ');

MLdecision = strings(numedges,1);
MAPdecision = strings(numedges,1);
pfalse_alarm = zeros(numedges,1);
pmiss = zeros(numedges,1);
likelihoodtables = cell(numedges,1);

for i = 1:numedges
    
    Xvals = (0:n)';
    PxH1 = binopdf(Xvals,n,p1);
    PxH0 = binopdf(Xvals,n,p0);
    T = table(Xvals, PxH1, PxH0, 'VariableNames', {'X','P_X_given_H1','P_X_given_H0'});
    likelihoodtables{i} = T;
    fprintf('\nLikelihood Table for Edge %s:\n', edgenames{i});
    disp(T);
    xi = X(i);
    
    if PxH1(xi+1) > PxH0(xi+1)
        MLdecision(i) = 'H1';
    else
        MLdecision(i) = 'H0';
    end
    
    dec_ML = PxH1 > PxH0;
    pfalse_alarm(i) = sum(PxH0(dec_ML));
    pmiss(i) = sum(PxH1(~dec_ML));
    
    postH1 = PxH1(xi+1)*prior_H1;
    postH0 = PxH0(xi+1)*(1-prior_H1);
    if postH1 > postH0
        MAPdecision(i) = 'H1';
    else
        MAPdecision(i) = 'H0';
    end
end

finalTable = table(X, MLdecision, MAPdecision, pfalse_alarm, pmiss, 'RowNames', edgenames, 'VariableNames', {'Successes','ML','MAP','P_false_alarm','P_miss'});

disp('Final ML and MAP Decisions Table:');
disp(finalTable);
fprintf('\nAnalysis complete. Edge-level ML and MAP decisions computed.\n');
