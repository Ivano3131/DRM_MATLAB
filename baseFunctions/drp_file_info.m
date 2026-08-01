function info = drp_file_info(file)
% drp_file_info  Read the header of a .drp file without loading the data.
%
%   info = drp_file_info("map.drp")
%
% A .drp is the single-file export of the DRM acquisition software.  It is a
% flat little-endian binary with no magic number and no text, laid out as
%
%   bytes  0 .. 15      header, 8 x uint16
%   bytes 16 .. 16+2P   preview image, P = width*height uint16, ROW major
%   the rest            the DRPs, PIXEL major: for every pixel in turn, all
%                       num_img uint16 intensities of that pixel's DRP
%
% and the header is
%
%   1 ph_num   2 th_num   3 ph_max   4 ph_min
%   5 th_max   6 th_min   7 width    8 height        (angles in whole degrees)
%
% The pixel-major body is the whole point of the format: a viewer can seek
% straight to one pixel's DRP instead of touching 1440 images.  It is also why
% loading is done a row at a time in drp_file_loader - an ROI is contiguous
% within a row and scattered across rows.
%
% Within a pixel the num_img values run in the SAME order as the phitheta
% table drp_loader builds from a folder, i.e. azimuth fastest:
%
%       k = th_index * ph_num + ph_index,   ph_index = phi / ph_step
%                                           th_index = (theta-th_min)/th_step
%
% (verified against the matching JPEG folder: 102_35.jpg lands on k = 634 =
% 5*120+34 with correlation 0.995, 180_40.jpg on k = 780 = 6*120+60.)
%
% OUTPUT ------------------------------------------------------------------
%   info.ph_num .. info.height   the eight header fields, as double
%   info.ph_step, info.th_step   derived degree steps
%   info.num_img, info.num_pix   ph_num*th_num, width*height
%   info.dataOffset              byte offset of the first pixel's DRP
%   info.previewOffset           byte offset of the preview image (16)
%   info.fileBytes               actual size on disk
%   info.expectedBytes           size the header implies
%   info.exp_para                ready-made struct of the six angle fields,
%                                so an experiment can be described BY the file
%                                rather than typed in beside it
%
% See also drp_file_loader, drp_loader.
% -------------------------------------------------------------------------
    arguments
        file (1,1) string
    end

    if ~isfile(file)
        error('drp_file_info:notFound','DRP file not found: %s', file);
    end

    % 'ieee-le' rather than the machine default: the format is little-endian
    % wherever it is read.
    fid = fopen(file, 'r', 'ieee-le');
    if fid < 0
        error('drp_file_info:cannotOpen','Cannot open DRP file: %s', file);
    end
    closeFid = onCleanup(@() fclose(fid));

    hdr = fread(fid, 8, '*uint16');
    if numel(hdr) < 8
        error('drp_file_info:shortHeader', ...
            '%s is only %d bytes - too short to hold a 16-byte DRP header.', ...
            file, numel(hdr)*2);
    end
    hdr = double(hdr);

    info.file   = file;
    info.ph_num = hdr(1);
    info.th_num = hdr(2);
    info.ph_max = hdr(3);
    info.ph_min = hdr(4);
    info.th_max = hdr(5);
    info.th_min = hdr(6);
    info.width  = hdr(7);
    info.height = hdr(8);

    if any([info.ph_num info.th_num info.width info.height] < 1)
        error('drp_file_info:badHeader', ...
            ['Header of %s does not look like a DRP header ' ...
             '(ph_num=%d th_num=%d width=%d height=%d).  ' ...
             'Is this really a .drp export?'], ...
            file, info.ph_num, info.th_num, info.width, info.height);
    end

    info.ph_step = local_step(info.ph_max, info.ph_min, info.ph_num);
    info.th_step = local_step(info.th_max, info.th_min, info.th_num);

    info.num_img = info.ph_num * info.th_num;
    info.num_pix = info.width  * info.height;

    info.previewOffset = 16;
    info.dataOffset    = 16 + info.num_pix * 2;
    info.expectedBytes = info.dataOffset + info.num_pix * info.num_img * 2;

    d = dir(file);
    info.fileBytes = d.bytes;

    % The six numbers an exp_para needs, straight from the file.
    info.exp_para = struct( ...
        'th_min', info.th_min, 'th_max', info.th_max, 'th_num', info.th_num, ...
        'ph_min', info.ph_min, 'ph_max', info.ph_max, 'ph_num', info.ph_num);
end

% =========================================================================
function s = local_step(vmax, vmin, n)
% One sample means no step; the file still describes a valid single slice.
if n > 1
    s = (vmax - vmin) / (n - 1);
else
    s = 0;
end
end
