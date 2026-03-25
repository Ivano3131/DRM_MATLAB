function xyz_Point = plotDataIPF(data,orientations,vecProj,options)
    % plot datapoints projecting on inverse pole figure
    % data contains the colorcoding for every datapoint
    % vecProj is the projecting direction of the orientations onto the IPF
    arguments
        data
        orientations
        vecProj
        %options.cs = crystalSymmetry('cubic')
        options.cs = {
          crystalSymmetry('6/mmm', [2.95, 2.95, 4.68], 'X||a*', 'Y||b', ...
          'mineral', 'Ti-Hex', 'color', 'light gray')
        };
        options.cmap = 'palura'
        options.linewidth (1,1) double = 2
        options.markersize (1,1) double = 8
        options.figSize (1,4) double = [2 2 5 5]
    end
    
    % find the outer countor projecting on fundemental region
    sR = fundamentalSector(options.cs);
    h = plotS2Grid(sR);
    xyz = abs(h.xyz);
    xy = sort(xyz(:,1:2),2,'descend');
    k = convhull(xy);
    boundaryX = xy(k,1);
    boundaryY = xy(k,2);
    
    % calculate the projection of different orientations
    r = inv(orientations) * vecProj;
    xyz_Point = normr(sort(abs(r.xyz),2,'ascend'));
    xy_Point = sort(xyz_Point(:,1:2),2,'descend');
    % figure("Units","centimeters","Position",options.figSize)
    scatter(xy_Point(:,1),xy_Point(:,2),options.markersize,data,'filled');
    set(gca,'Visible','off')
    axis equal
    hold on
    plot(boundaryX,boundaryY,'k',LineWidth=options.linewidth)
    hold off

end