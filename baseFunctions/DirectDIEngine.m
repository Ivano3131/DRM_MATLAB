function indexResult = DirectDIEngine(drpM, drpLib, options)
    arguments
        drpM % measured
        drpLib % library
        options.K (1,1) double = 1
    end
    [n1,n2] = size(drpM);
    EUmap = zeros(3,n1,n2);
    drplist_s = zeros(length(drpLib.drpDic), numel(drpM{1,1}));
    for ii = 1:length(drpLib.drpDic)
        drplist_s(ii,:) = double(reshape(drpLib.drpDic{ii},1,[]))/256;
    end
    for ii = 1: n1
        drplist_m = zeros(n2,numel(drpM{1,1}));
        for jj = 1:n2
            drplist_m(jj,:) = double(reshape(drpM{ii,jj},1,[]))/256;
        end

        [Idx, D] = knnsearch(drplist_s, drplist_m, K=options.K); %search the full length of the DRPs

        EUmap(:,ii,:) = drpLib.eulerDic(Idx(:,1),:)';
        indexResult.idxMap(ii,:) = Idx(:,1);
        workbar(ii/n1, sprintf("processing %d / %d lines...",[ii n1]));
    end
    indexResult.euMap = permute(EUmap, [2,3,1]);
end