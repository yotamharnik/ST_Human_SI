function [barcodes] = read_spots2exclude(path2file)
% this function gets a csv output from lupa and takes only the barcodes
t = readtable(path2file);
barcodes = table2cell(t(:,1));
for i=1:length(barcodes)
        str=barcodes{i};
        str(findstr(str,'-'))='_';
        barcodes{i}=str;
end
end