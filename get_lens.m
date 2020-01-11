function [er,mur] = get_lens(WG0,l_rwg,neff,er,mur,R,er0,mur0,isPEC,b1,b2,b21,b22,flag)
% hemispherical lense is hard code - In addition it is not a good lens
a2 = b2 - b1;
RL = 2*20;
rho = sqrt(R^2+a2^2/4);
phi = asin(a2/(2*rho));
theta = linspace(-phi,phi,200);
face0 = linspace(b21,b22,200);
if strcmp(flag,'hem')
    shift = 0;
    for s = 1:40
        if s < 35
            shift = shift+1;
        end
        er(280 +s,119-shift:121+shift) = neff*er0;
        mur(280 +s,119-shift:121+shift) = neff*mur0;
    end
end
if (strcmp(flag,'grad') || strcmp(flag,'lens'))
    % Start getting the graded index
    for s = 1:length(theta)
        [n{s},pos(s,:,:)] = get_ref(WG0,l_rwg,R,theta(s),rho,RL,face0(s)); % Required Refractive index needed in every PHYSICAL point
    end
    
    for s = 1:length(theta)
        cell_neff = pos2cell(pos(s,:,:));
        % When non-unique cells appear take avg from corresponding neff    
        [~, ~, ib] = unique(cell_neff, 'rows');
        numoccurences = accumarray(ib, 1);
        indices = accumarray(ib, find(ib), [], @(rows){rows});
        % from cell_neff assign dielectric and magnetic constants
        for q = 1:length(indices)
            ind = indices{q}(1);
            er(cell_neff(ind,1),cell_neff(ind,2)) = er0*mean(n{s}(indices{q})).^2;
%             mur(cell_neff(ind,1),cell_neff(ind,2)) = mur0*mean(n{s}(indices{q})); % One could have permitivity and permiability to get better matching
        end
    end
end

