function cell = pos2cell(pos)
global dx
for s = 1:length(pos)
    cell(s,:) = pos(:,:,s)./dx;
end
cell = round(cell);

