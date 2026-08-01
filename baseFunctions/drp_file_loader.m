function [igrey,phitheta,pos,img_sample] = drp_file_loader(exp_para,pos,options)
% drp_file_loader  Load a single .drp export into the same stack drp_loader
% returns from a folder of images.
%
%   [igrey,phitheta,pos,img_sample] = drp_file_loader(exp_para,pos,file="map.drp")
%
% Drop-in for the folder path: identical outputs, so igray2drp and everything
% downstream is unchanged.
%
%   igrey       n1 x n2 x num_img uint8 (or uint16, see precision)
%   phitheta    num_img x 2 [phi theta] in degrees, in stack order
%   pos         the ROI actually used
%   img_sample  the cropped preview image, for overlays and clicking
%
% See drp_file_info for the byte layout.  Two things differ from the folder
% path and both follow from the file already being a finished export:
%
%   * there is no `scale`.  A folder holds full-resolution camera frames that
%     you shrink on the way in; a .drp was written at whatever resolution the
%     acquisition software exported, and `pos` is in THOSE pixels.  For the
%     Iven cropped_left pair the .drp is 859 x 963, which is exactly what the
%     folder path produces at scale 0.4458 - so pos carries over unchanged.
%   * there is no filename parsing, so no dependence on the <phi>_<theta>
%     naming convention and no way for a missing or misnamed image to shift
%     the stack silently.  The angles come from the header.
%
% OPTIONS -----------------------------------------------------------------
%   file        path to the .drp.  "" opens a file chooser.
%   precision   "uint8" (default) divides the stored uint16 by 257, giving
%               exactly the dynamic range the folder path gives, so existing
%               thresholds (the index_num > 3e4 background mask) still hold.
%               "uint16" keeps the full depth the file actually carries -
%               ~6900 distinct levels in the Ti64 data rather than 256.
%               Indexing normalises every DRP anyway, so this changes
%               quantisation noise, not the scale of the score.
%   strict      true (default) errors if the header disagrees with exp_para.
%               The dictionary is built from exp_para BEFORE the data is
%               loaded, so a disagreement means the two are not comparable -
%               see drp_file_info(file).exp_para to set exp_para from the file.
%   verbose     true (default) prints the geometry and a progress bar.
%
% DOES IT GIVE THE SAME MAP AS THE FOLDER? --------------------------------
% The two paths load the same measurement, but they are NOT the same numbers,
% and on this dataset they do not produce the same map.  Measured on the
% cropped_left pair (80x80 ROI, 5 deg dictionary):
%
%   spatial alignment      identical - the frame correlation peaks at shift
%                          (0,0) at 0.989 and falls to 0.91 one pixel away
%   intensity              linear, drp/257 = 0.886*jpg + 1.9, rms 4.2 levels
%   assembled DRPs         median correlation 0.98 between the two paths
%   fit to the cube model  median best-match score 0.856 folder, 0.848 .drp
%   SAME dictionary entry  only 35% of pixels, misorientation median 11 deg
%
% The last line is not a loading error - it is the model's own ambiguity.  A
% typical pixel has nine dictionary candidates within 1% of its best score
% (best 0.8557, second 0.8554), so the ~2% by which JPEG compression and
% MATLAB's resize kernel differ from the raw export is enough to reorder them.
% Jittering the folder stack by +-1 grey level, a perturbation ten times
% smaller, leaves 99.3% of pixels on the same entry - the match is stable, the
% top of the dictionary is just crowded.  Trust res.score and the ambiguity
% report from cubeAmbiguity, and do not expect the two inputs to be
% interchangeable pixel by pixel.
%
% MEMORY ------------------------------------------------------------------
% The stack is n1*n2*num_img bytes at uint8, twice that at uint16 - 1.2 GB for
% a full 859 x 963 x 1440 map.  The file is read one ROI row at a time and
% converted immediately, so peak memory is the stack plus one row, not the
% whole 2.4 GB file.  Crop with `pos` if that is too much.
%
% See also drp_file_info, drp_loader, igray2drp.
% -------------------------------------------------------------------------
    arguments
        exp_para struct
        pos (1,4) double
        options.file (1,:) string = ""
        options.precision (1,1) string ...
            {mustBeMember(options.precision,["uint8","uint16"])} = "uint8"
        options.strict (1,1) logical = true
        options.verbose (1,1) logical = true
    end

    file = options.file;
    if file == ""
        [f,p] = uigetfile({'*.drp','DRP file (*.drp)'},'Select a .drp file');
        if isequal(f,0)
            error('drp_file_loader:cancelled','No .drp file selected.');
        end
        file = string(fullfile(p,f));
    end

    info = drp_file_info(file);

    % ---- does the file describe the experiment we built a dictionary for? -
    local_checkGeometry(info, exp_para, options.strict, options.verbose);

    % ---- is the file complete? -------------------------------------------
    % Short means truncated, and a truncated pixel-major file is not simply a
    % smaller map: the missing bytes are the LAST pixels, so the rows that do
    % load are still correct.  Say how far it got rather than guessing.
    if info.fileBytes < info.expectedBytes
        gotPix = floor((info.fileBytes - info.dataOffset) / (info.num_img*2));
        error('drp_file_loader:truncated', ...
            ['%s is truncated: %d bytes on disk, %d expected for a ' ...
             '%d x %d map of %d images.  Only %d of %d pixel records are ' ...
             'present (%.1f%%).'], ...
            file, info.fileBytes, info.expectedBytes, ...
            info.width, info.height, info.num_img, ...
            max(gotPix,0), info.num_pix, 100*max(gotPix,0)/info.num_pix);
    elseif info.fileBytes > info.expectedBytes && options.verbose
        fprintf(['drp_file_loader: %d trailing bytes after the last DRP, ' ...
                 'ignored.\n'], info.fileBytes - info.expectedBytes);
    end

    if options.verbose
        fprintf('drp_file_loader: %s\n', file);
        fprintf(['  %d x %d pixels, %d images ' ...
                 '(phi %g:%g:%g, theta %g:%g:%g)\n'], ...
            info.width, info.height, info.num_img, ...
            info.ph_min, info.ph_step, info.ph_max, ...
            info.th_min, info.th_step, info.th_max);
    end

    fid = fopen(file, 'r', 'ieee-le');
    if fid < 0
        error('drp_file_loader:cannotOpen','Cannot open DRP file: %s', file);
    end
    closeFid = onCleanup(@() fclose(fid));

    % ---- preview image ----------------------------------------------------
    % Stored row major, so fill [width height] and transpose - MATLAB's
    % reshape runs down columns.
    fseek(fid, info.previewOffset, 'bof');
    preview = fread(fid, info.num_pix, '*uint16');
    preview = reshape(preview, info.width, info.height).';

    % ---- ROI --------------------------------------------------------------
    if isequal(pos, zeros(1,4))
        figure(101)
        imshow(im2uint8(preview) * 2, 'Border', 'tight')
        roi = drawrectangle;
        pos = roi.Position;
        close 101
    end

    % Rather than re-deriving imcrop's rounding rules, crop two coordinate
    % grids with imcrop itself.  Whatever rectangle the folder path would have
    % produced, this produces the same one - including the clamping when pos
    % runs past the edge.
    rGrid = repmat((1:info.height)', 1, info.width);
    cGrid = repmat(1:info.width, info.height, 1);
    rows  = imcrop(rGrid, pos);
    cols  = imcrop(cGrid, pos);
    if isempty(rows)
        error('drp_file_loader:emptyROI', ...
            'ROI [%g %g %g %g] selects no pixels of the %d x %d map.', ...
            pos(1), pos(2), pos(3), pos(4), info.width, info.height);
    end
    rows = rows(:,1).';
    cols = cols(1,:);
    n1   = numel(rows);
    n2   = numel(cols);

    % The single-fread-per-row scheme below assumes the ROI columns are one
    % consecutive run.  imcrop always gives that; check it rather than trust it,
    % because a gap would silently shift every pixel after it.
    if n2 > 1 && ~isequal(cols, cols(1):cols(end))
        error('drp_file_loader:nonContiguousROI', ...
            'ROI columns are not consecutive - cannot read rows in one block.');
    end

    img_sample = im2uint8(preview(rows, cols));

    if options.verbose
        fprintf('  ROI [%g %g %g %g] -> %d x %d pixels, %s stack (%.2f GB)\n', ...
            pos(1), pos(2), pos(3), pos(4), n1, n2, options.precision, ...
            n1*n2*info.num_img*(1 + (options.precision=="uint16"))/2^30);
    end

    % ---- phitheta, in stack order ----------------------------------------
    % Built exactly as drp_loader builds it, because the file's within-pixel
    % order is exactly that order.  splitDRP reads this table to place each
    % sample, so the two paths produce identical DRPs.
    phitheta = zeros(info.num_img, 2);
    for ii = 1:info.num_img
        phitheta(ii,1) = rem(ii-1, info.ph_num) * info.ph_step;
        phitheta(ii,2) = floor((ii-1) / info.ph_num) * info.th_step + info.th_min;
    end

    % ---- the stack --------------------------------------------------------
    % One fread per ROI row: within a row the wanted pixels are adjacent, and
    % each pixel's num_img values are contiguous, so cols(1)..cols(end) is a
    % single run of n2*num_img uint16.
    igrey  = zeros(n1, n2, info.num_img, options.precision);
    rowLen = info.width * info.num_img * 2;     % bytes per image row
    base   = info.dataOffset + (cols(1)-1) * info.num_img * 2;

    for ii = 1:n1
        fseek(fid, base + (rows(ii)-1) * rowLen, 'bof');
        blk = fread(fid, [info.num_img, n2], '*uint16');
        if size(blk,2) < n2
            error('drp_file_loader:readFailed', ...
                'Short read on image row %d of %s.', rows(ii), file);
        end
        if options.precision == "uint8"
            blk = im2uint8(blk);                % /257, full-range
        end
        igrey(ii,:,:) = reshape(blk.', 1, n2, info.num_img);
        if options.verbose
            workbar(ii/n1, sprintf('Reading DRPs, %d / %d rows',[ii n1]));
        end
    end
end

% =========================================================================
function local_checkGeometry(info, exp_para, strict, verbose)
% exp_para drives the dictionary and splitDRP's binning; the file drives the
% data.  Where they disagree the map would be silently wrong, so say so here.
fields = ["th_min","th_max","th_num","ph_min","ph_max","ph_num"];
bad    = strings(0,1);
for f = fields
    if isfield(exp_para, f) && exp_para.(f) ~= info.(f)
        bad(end+1) = sprintf('  %-7s exp_para %g, file %g', ...
            f, exp_para.(f), info.(f)); %#ok<AGROW>
    end
end
if isempty(bad)
    return
end
msg = sprintf(['The .drp header does not match exp_para:\n%s\n' ...
    'The dictionary was built from exp_para, so the two are not ' ...
    'comparable.  Set exp_para from the file with\n' ...
    '    hdr = drp_file_info(file);  exp_para = ... hdr.exp_para fields\n' ...
    'or pass strict=false to load anyway.'], strjoin(bad, newline));
if strict
    error('drp_file_loader:geometryMismatch', '%s', msg);
elseif verbose
    warning('drp_file_loader:geometryMismatch', '%s', msg);
end
end
