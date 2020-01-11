function [phase,freq] = cst_txt_ff(fid)
no = 0;
while ~feof(fid)
    tline = fgetl(fid);
    if contains(tline,'Theta')
        no = no+1;
        count = 1;
    end    
    tline = str2num(tline);
    if ~isempty(tline)
        freq(no,count) = tline(3);
        phase(no,count) = tline(2);
        count = count + 1;
    end
end
% freq = freq./1e9;