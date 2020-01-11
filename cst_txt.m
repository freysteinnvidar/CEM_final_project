function [phase,freq] = cst_txt(fid)
no = 0;
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline,'e-field')
        no = no+1;
        count = 1;
    end    
    tline = str2num(tline);
    if ~isempty(tline)
        freq(no,count) = tline(2);
        phase(no,count) = tline(1);
        count = count + 1;
    end
end
% freq = freq./1e9;