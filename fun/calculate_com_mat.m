function [com_vec] = calculate_com_mat(mat,com_lim,REMOVE_MIN)
if nargin < 2
    com_lim = 1;
end
if nargin<3
    REMOVE_MIN=0;
end
com_vec = zeros(size(mat,1),1);
for i = 1 :size(mat,1)
    vec=(mat(i,:));
    if REMOVE_MIN
        vec=(vec-min(vec))/(max(vec)-min(vec));
    end
    com_vec(i) = nansum(linspace(0,com_lim,length(vec)).*vec)/nansum(vec);
end
end