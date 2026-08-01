function export_ctf(ctfFile, EUmap, mask, opts)
% Write a predicted orientation map as an Oxford/HKL Channel Text File (.ctf)
%
% The MATLAB twin of the python export_ctf: same header, same row-major data
% block (X fastest, then Y), so the file drops straight into MTEX, Channel or
% a ctf_compare script.
%
% INPUT -------------------------------------------------------------------
%   ctfFile   output path, e.g. "pred.ctf"
%   EUmap     n1 x n2 x 3 [phi1 Phi phi2] IN DEGREES, as returned by
%             matchDRPcube / IndexingEngine.  Row = Y, column = X.
%   mask      n1 x n2 logical/numeric sample mask; 0 = not indexed.  Those
%             pixels are written as Phase 0 with zero Euler angles, exactly
%             like the python version.  Default: every pixel indexed.
%
% OPTIONS -----------------------------------------------------------------
%   xStep,yStep  step size in um (default 1).
%   phase        struct with fields abc, ang, name, laue, sg.  Default is
%                alpha-Ti (hcp, Laue 9, spacegroup 194).  For prior-beta use
%                struct('abc',[3.32 3.32 3.32],'ang',[90 90 90], ...
%                       'name','Ti-BCC','laue',11,'sg',229) - writing the
%                wrong one makes every reader apply the wrong symmetry.
%   author       Author line (default "MATLAB DRM export").
% -------------------------------------------------------------------------
arguments
    ctfFile   (1,1) string
    EUmap     double
    mask                = []
    opts.xStep (1,1) double = 1
    opts.yStep (1,1) double = 1
    opts.phase struct   = struct('abc',[2.95 2.95 4.68],'ang',[90 90 120], ...
                                 'name','Ti-Hex','laue',9,'sg',194)
    opts.author (1,1) string = "MATLAB DRM export"
end

[nY, nX, nA] = size(EUmap);
if nA ~= 3, error('EUmap must be n1 x n2 x 3 (phi1, Phi, phi2).'); end
if isempty(mask), mask = true(nY,nX); end
if ~isequal(size(mask),[nY nX])
    error('Mismatch between map size (%d x %d) and mask size (%d x %d).', ...
        nY, nX, size(mask,1), size(mask,2));
end

rowmaj = @(A) reshape(double(A).',[],1);      % matrix -> CTF row-major vector

phase = rowmaj(mask ~= 0);
eul   = [rowmaj(EUmap(:,:,1)), rowmaj(EUmap(:,:,2)), rowmaj(EUmap(:,:,3))];
eul(~isfinite(eul)) = 0;
eul(phase == 0, :)  = 0;                      % zero out missing points

[xg, yg] = meshgrid((0:nX-1)*opts.xStep, (0:nY-1)*opts.yStep);
x = rowmaj(xg);  y = rowmaj(yg);

fid = fopen(ctfFile,'wt');
if fid < 0, error('Cannot open %s for writing.',ctfFile); end
cl = onCleanup(@() fclose(fid)); %#ok<NASGU>

[~,prj] = fileparts(char(ctfFile));
fprintf(fid,'Channel Text File\n');
fprintf(fid,'Prj\t%s\n',prj);
fprintf(fid,'Author\t%s\n',opts.author);
fprintf(fid,'JobMode\tGrid\n');
fprintf(fid,'XCells\t%d\n',nX);
fprintf(fid,'YCells\t%d\n',nY);
fprintf(fid,'XStep\t%.4f\n',opts.xStep);
fprintf(fid,'YStep\t%.4f\n',opts.yStep);
fprintf(fid,'AcqE1\t0\nAcqE2\t0\nAcqE3\t0\n');
fprintf(fid,['Euler angles refer to Sample Coordinate system (CS0)!\t' ...
             'Mag\t100\tCoverage\t100\tDevice\t0\tKV\t20\t' ...
             'TiltAngle\t70\tTiltAxis\t0\n']);
fprintf(fid,'Phases\t%d\n',numel(opts.phase));
for k = 1:numel(opts.phase)
    p = opts.phase(k);
    fprintf(fid,'%.3f;%.3f;%.3f\t%.3f;%.3f;%.3f\t%s\t%d\t%d\n', ...
        p.abc(1),p.abc(2),p.abc(3), p.ang(1),p.ang(2),p.ang(3), ...
        p.name, p.laue, p.sg);
end
fprintf(fid,'Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n');

M = [phase, x, y, eul].';
fprintf(fid,'%d\t%.2f\t%.2f\t0\t0\t%.4f\t%.4f\t%.4f\t0\t0\t0\n',M);

fprintf('CTF export complete: %s  (%d x %d px, %d indexed)\n', ...
    ctfFile, nX, nY, sum(phase ~= 0));
end
