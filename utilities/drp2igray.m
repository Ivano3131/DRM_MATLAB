function igray = drp2igray(drp_original,exp_para)
% convert cell-formated DRP back to stack of raw images
    arguments
        drp_original cell
        exp_para struct
    end
    [n1,n2] = size(drp_original);
    % 'like' the DRPs themselves: a stack loaded from a .drp at precision
    % "uint16" would otherwise saturate on the way back out.
    igray = zeros(n1,n2,exp_para.ph_num*exp_para.th_num,'like',drp_original{1,1});
    for ii = 1:n1
        for jj = 1:n2
            drp_temp = drp_original{ii,jj};
            igray(ii,jj,:) = reshape(drp_temp,1,[]);
        end
    end
    fprintf("transformation from cell to matrix finished!\n")
end