function [S_selected] = ratio_selection(scores,ratio,type)
% Select top specific ratio of data

% scores = Z_sample_normal;
% ratio = 0.1;

N = size(scores,1);
NR = randperm(N);
S = [scores,(1:1:N)'];
SR = S(NR,:);
if type == 1
    S_sort = sortrows(SR,1,'ascend');
else
    S_sort = sortrows(SR,1,'descend');
end
S_selected = S_sort(1:round(N*ratio),:);





