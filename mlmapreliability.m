clc; clear;

N = input('Enter number of nodes: ');
E = input('Enter number of edges: ');
fprintf('Enter edges as pairs [u v] (one per row):\n');
edges = zeros(E,2);
for i = 1:E
    edges(i,:) = input(sprintf('Edge %d: ', i));
end

edgeNames = strcat(string(edges(:,1)),'-',string(edges(:,2)));
fprintf('Enter ML estimated probabilities for each edge:\n');
p_ML = zeros(E,1);
for i = 1:E
    p_ML(i) = input(sprintf('Edge %s ML probability: ', edgeNames{i}));
end

fprintf('Enter MAP estimated probabilities for each edge:\n');
p_MAP = zeros(E,1);
for i = 1:E
    p_MAP(i) = input(sprintf('Edge %s MAP probability: ', edgeNames{i}));
end

s = 1;       
t = N;       
cutSets = {}; 
allEdgesIdx = 1:E;

for k = 1:E
    subsets = nchoosek(allEdgesIdx,k);
    for i = 1:size(subsets,1)
        tempEdges = edges;
        tempEdges(subsets(i,:),:) = [];
        G = graph(tempEdges(:,1), tempEdges(:,2));
        if ~ismember(t, bfsearch(G,s))
            cutSets{end+1} = sort(subsets(i,:)); 
        end
    end
end


cutSets = unique(cellfun(@(x) num2str(x), cutSets, 'UniformOutput', false));
cutSets = cellfun(@(x) str2num(x), cutSets, 'UniformOutput', false); 

fprintf('\nFound %d cut sets separating node %d and %d\n', length(cutSets), s, t);

R_ML = 1; R_MAP = 1;
failure_ML = 0; failure_MAP = 0;

for k = 1:length(cutSets)
    idx = cutSets{k};
    pf_ML = prod(1 - p_ML(idx));
    pf_MAP = prod(1 - p_MAP(idx));
    failure_ML = failure_ML + pf_ML;
    failure_MAP = failure_MAP + pf_MAP;
end

R_ML = max(0, 1 - failure_ML);
R_MAP = max(0, 1 - failure_MAP);

fprintf('\nUnion-bound estimate of reliability (node %d -> node %d):\n', s, t);
fprintf('ML-based estimate: %.4f\n', R_ML);
fprintf('MAP-based estimate: %.4f\n', R_MAP);
