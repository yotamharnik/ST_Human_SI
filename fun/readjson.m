function val = readjson(fname)
%fname = 'Z:\yotamhar\NGS\211229_A00929_0540_BHVK7NDRXY_Visium_H-duo_P1-P9\BHVK7NDRXY\P1\outs\spatial\scalefactors_json.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);