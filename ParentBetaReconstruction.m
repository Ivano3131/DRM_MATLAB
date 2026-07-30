%% Parent beta grain reconstruction from alpha EBSD (Ti-6Al-4V, Widmanstatten)
%  Reconstructs the prior-beta grain structure from a child alpha map using the
%  Burgers OR and MTEX's variant-graph algorithm, plots the beta grain map, and
%  exports the reconstructed beta orientations as a .ctf file on the SAME grid
%  as the input map (so it stays pixel-aligned with the alpha map / DRM).
%
%  Requires MTEX 5.8 or newer (for calcVariantGraph / clusterVariantGraph).
%  Things you must adapt to your data are marked  <-- ADAPT.

%% 1. Load your EBSD data ------------------------------------------------------
% Replace this with your own import line. Keep the reference-frame conversion
% consistent with how you imported the map you co-registered against the DRM.

%addpath("C:\Users\mrbla\Desktop\mrbla-downloads\Relevant Downloads\MTEX Code\mtex-6.1.0")
%startup_mtex
%ebsd = EBSD.load("C:\Users\mrbla\Desktop\Cambridge\Thesis\Results\Ti64-BottomRightScratch.ctf",'convertEuler2SpatialReferenceFrame');   % <-- ADAPT
inFile = "C:\Users\mrbla\Desktop\Cambridge\Thesis\Results\AM-Ti64\EBSD-AMTi64.ctf";
ebsd   = EBSD.load(inFile, 'convertEuler2SpatialReferenceFrame');

% >>> NEW: pristine, full-grid copy taken BEFORE any pixels are deleted.
%     The export in section 11 uses this so the output .ctf has exactly the same
%     XCells / YCells / XStep / YStep as the input file.
ebsd0 = ebsd;

% Name of the alpha (child) phase exactly as it appears in the ebsd printout:
alphaName = 'Ti-Hex';   % <-- ADAPT (e.g. 'Ti-Hex', 'Titanium', 'Ti (Alpha)')

%% 2. Define the beta (parent) phase -------------------------------------------
% Only the cubic *symmetry* matters for the reconstruction; the lattice constant
% affects direction plots only, so its exact value is not critical.
betaName = 'Titanium cubic';

%% 3. Reconstruct alpha grains
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',5*degree); %3 for wrought

% OPTIONAL cleanup — remove the PIXELS of tiny grains, then recompute.
% (100 px is very aggressive on fine Widmanstatten laths; keep it small.)
ebsd(grains(grains.numPixel < 3)) = []; % 10 for wrought
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'threshold',5*degree); %3 for wrought

%% 4. Set up the reconstructor
csBeta = crystalSymmetry('m-3m',[3.32 3.32 3.32],'mineral',betaName);
job = parentGrainReconstructor(ebsd,grains, csBeta);
job.p2c = orientation.Burgers(job.csParent,job.csChild);   % defines the beta parent

% (optional) refine the OR to the data. For clean Widmanstatten Ti64 this is
% often unnecessary; uncomment if the alpha-alpha boundaries deviate from BOR.
job.calcParent2Child;

%% 5. Variant-graph reconstruction --------------------------------------------
job.calcVariantGraph('threshold',1.5*degree, 'tolerance', 1.5*degree);   % build the graph
job.clusterVariantGraph('numIter',3);            % cluster -> parent probabilities
job.calcParentFromVote('minFit', 2.5*degree);    % assign most-likely parent

%% 6. Clean up -----------------------------------------------------------------
job.mergeSimilar('threshold',5*degree);   % merge fragments of the same beta grain
job.mergeInclusions('maxSize',10);        % absorb tiny unassigned inclusions

%% 7. Plot the parent beta grain map ------------------------------------------
ipfKey = ipfColorKey(job.csParent);
ipfKey.inversePoleFigureDirection = vector3d.Z;   % map-normal IPF; change if needed
figure; plot(ipfKey);

color = ipfKey.orientation2color(job.parentGrains.meanOrientation);
figure; plot(job.parentGrains,color,'figSize','large');
title('Reconstructed prior-\beta grains');

%% 8. Overlay reconstructed beta boundaries on the alpha map -------------------
figure;
plot(ebsd(job.csChild),ebsd(job.csChild).orientations,'figSize','large');
hold on
betaGrains = smooth(job.parentGrains,5);
plot(betaGrains.boundary,'lineWidth',3,'lineColor','white');
hold off
title('Alpha IPF with reconstructed prior-\beta boundaries');

%% 9. Per-pixel reconstructed beta orientations --------------------------------
% Newer MTEX exposes this as job.ebsd; older versions used job.calcParentEBSD.
parentEBSD = job.ebsd;
figure;
plot(parentEBSD(job.csParent), ...
    ipfKey.orientation2color(parentEBSD(job.csParent).orientations),'figSize','large');
title('Per-pixel reconstructed \beta orientation');

%% 10. Test the point-3 hypothesis --------------------------------------------
% Do the flat and rough regions reconstruct to DIFFERENT beta parents, despite
% their alpha orientations differing by only ~4 deg?
pos_flat  = [ 10 ,  50 ];    % <-- ADAPT
pos_rough = [ 6,  27];       % <-- ADAPT
pg = job.parentGrains;

id_flat  = pg.findByLocation(pos_flat);
id_rough = pg.findByLocation(pos_rough);

pick = @(id) selectGrain(pg,id);
g_flat  = pick(id_flat);
g_rough = pick(id_rough);

fprintf('flat:  grain id %d, area %.1f um^2\n', g_flat.id,  g_flat.area);
fprintf('rough: grain id %d, area %.1f um^2\n', g_rough.id, g_rough.area);
fprintf('beta parent misorientation: %.1f deg\n', ...
    angle(g_flat.meanOrientation, g_rough.meanOrientation)/degree);

% Compare this to the ~4 deg between the alpha orientations. A large beta
% misorientation with small alpha misorientation is the signature you want.


%% =====================================================================
%% 11. EXPORT THE RECONSTRUCTED BETA MAP AS .ctf  (new)
%% =====================================================================
%  The output is written on the FULL original grid, in the CTF row-major order
%  (X fastest, then Y), so row k of the data block corresponds to exactly the
%  same pixel as row k of the input file. Pixels with no assigned beta parent
%  (unindexed, deleted small grains, unreconstructed) are written as Phase 0.

ctfOut = "C:\Users\mrbla\Desktop\Cambridge\Thesis\Results\AM-Ti64\EBSD-AMTi64-beta.ctf";  % <-- ADAPT

useGrainMean  = true;  % false -> per-pixel reconstructed beta orientation
                        % true  -> every pixel gets its parent grain's MEAN orientation
                        %          (cleaner ground truth for ML, no variant scatter)
undoImportRot = true;   % write Euler angles back in the input file's frame, i.e.
                        % undo 'convertEuler2SpatialReferenceFrame'. Section 12
                        % verifies this by reloading the file — check the printout.
carryOverBCBS = true;   % copy Bands/Error/MAD/BC/BS from the original map
                        % (handy as extra ML channels / quality masks)

assert(exist('ebsd0','var')==1, ...
    'ebsd0 is missing — add "ebsd0 = ebsd;" right after EBSD.load in section 1.');

% --- grid geometry, taken from the ORIGINAL map -------------------------------
eg = gridify(ebsd0);
[nY,nX] = size(eg);
assert(nX>1 && nY>1,'Map is not a 2-D grid.');
[Xg,Yg] = xyOf(eg);                    % nY x nX coordinate matrices
dx = Xg(1,2)-Xg(1,1);
dy = Yg(2,1)-Yg(1,1);
x0 = Xg(1,1);  y0 = Yg(1,1);
N  = nX*nY;

rowmaj = @(A) reshape(double(A).',[],1);   % matrix -> CTF row-major column vector

% --- pixels that carry a reconstructed beta orientation -----------------------
par  = job.ebsd(job.csParent);
oriP = par.orientations;

if useGrainMean
    try
        [tf,loc] = ismember(par.grainId, job.parentGrains.id);
        oriP(tf) = job.parentGrains.meanOrientation(loc(tf));
    catch ME
        warning('Grain-mean substitution failed (%s) — using per-pixel orientations.', ...
            ME.message);
    end
end

[px,py] = xyOf(par);
ix = round((px - x0)/dx) + 1;
iy = round((py - y0)/dy) + 1;
ok = ix>=1 & ix<=nX & iy>=1 & iy<=nY;
lin  = (iy(ok)-1)*nX + ix(ok);      % row-major linear index into the full grid
oriP = oriP(ok);

% --- Euler angles (Bunge, degrees) --------------------------------------------
oriW = oriP;
if undoImportRot
    oriW = rotation.byAxisAngle(xvector,180*degree) * oriW;
end
[e1,e2,e3] = Euler(oriW,'Bunge');
e1 = mod(e1/degree,360);  e2 = e2/degree;  e3 = mod(e3/degree,360);

% --- assemble the full-grid data block -----------------------------------------
d.phase = zeros(N,1);
d.x     = rowmaj(Xg);
d.y     = rowmaj(Yg);
d.bands = zeros(N,1);
d.err   = 3*ones(N,1);      % 3 = "not indexed" style error code
d.e1    = zeros(N,1);
d.e2    = zeros(N,1);
d.e3    = zeros(N,1);
d.mad   = zeros(N,1);
d.bc    = zeros(N,1);
d.bs    = zeros(N,1);

if carryOverBCBS
    d.bands = getProp(eg,'bands',rowmaj,N);
    d.mad   = getProp(eg,'mad',  rowmaj,N);
    d.bc    = getProp(eg,'bc',   rowmaj,N);
    d.bs    = getProp(eg,'bs',   rowmaj,N);
end

d.phase(lin) = 1;
d.err(lin)   = 0;
d.e1(lin)    = e1;
d.e2(lin)    = e2;
d.e3(lin)    = e3;

d.bands = round(d.bands);  d.bc = round(d.bc);  d.bs = round(d.bs);

% --- phase table ---------------------------------------------------------------
% To ALSO keep the alpha phase in the file, append a second entry here
% (3.232;3.232;5.147 / 90;90;120 / Ti-Hex / Laue 9 / SG 194) and set
% d.phase(alphaPixels) = 2 with the corresponding alpha Euler angles.
phases(1).abc  = [3.32 3.32 3.32];
phases(1).ang  = [90 90 90];
phases(1).name = betaName;
phases(1).laue = 11;    % m-3m
phases(1).sg   = 229;   % Im-3m (beta-Ti)

writeCTF(ctfOut,nX,nY,dx,dy,phases,d);

fprintf('\nWrote %s\n',ctfOut);
fprintf('  grid      : %d x %d px, step %.4f x %.4f um\n',nX,nY,dx,dy);
fprintf('  beta pixels: %d / %d  (%.1f %% coverage)\n', ...
    numel(lin),N,100*numel(lin)/N);

%% 12. Round-trip check ---------------------------------------------------------
% Reload the file the same way you loaded the input and compare orientations.
% max misorientation ~0 deg  -> undoImportRot is set correctly.
% max misorientation large   -> flip undoImportRot and re-export.
try
    chk  = EBSD.load(ctfOut,'convertEuler2SpatialReferenceFrame');
    chkP = chk(betaName);
    [cx,cy] = xyOf(chkP);
    lc = (round((cy-y0)/dy))*nX + round((cx-x0)/dx) + 1;
    [~,ia,ib] = intersect(lc,lin);
    dd = angle(chkP.orientations(ia),oriP(ib))/degree;
    fprintf('  round-trip: max %.4f deg, median %.4f deg over %d px\n', ...
        max(dd),median(dd),numel(dd));
    if max(dd) > 0.5
        warning(['Round-trip misorientation is large. Set undoImportRot = %d ' ...
                 'and re-run section 11.'],~undoImportRot);
    end
catch ME
    warning('Round-trip check skipped: %s',ME.message);
end


%% =====================================================================
%% local functions
%% =====================================================================
function g = selectGrain(pg,id)
    if isempty(id) || any(id <= 0)
        error('Position is not inside any reconstructed beta grain.');
    end
    g = pg(ismember(pg.id,id));      % id-based lookup
    if isempty(g), g = pg(id); end   % fallback if it returned a linear index
    g = g(1);
end

function [xx,yy] = xyOf(e)
% Coordinates, working across MTEX 5.x (e.x/e.y) and 6.x (e.pos)
    try
        xx = e.x;      yy = e.y;
    catch
        xx = e.pos.x;  yy = e.pos.y;
    end
end

function v = getProp(eg,name,rowmaj,N)
% Pull a gridded property, replacing gridify's NaN padding with 0
    if isfield(eg.prop,name)
        v = rowmaj(eg.prop.(name));
        v(~isfinite(v)) = 0;
    else
        v = zeros(N,1);
    end
end

function writeCTF(fname,nX,nY,dx,dy,phases,d)
% Minimal Channel Text File writer (Oxford/HKL .ctf), single or multi phase.
    fid = fopen(fname,'wt');
    if fid < 0, error('Cannot open %s for writing.',fname); end
    cl = onCleanup(@() fclose(fid)); %#ok<NASGU>

    [~,prj] = fileparts(char(fname));

    fprintf(fid,'Channel Text File\n');
    fprintf(fid,'Prj\t%s\n',prj);
    fprintf(fid,'Author\t[MTEX parentGrainReconstructor]\n');
    fprintf(fid,'JobMode\tGrid\n');
    fprintf(fid,'XCells\t%d\n',nX);
    fprintf(fid,'YCells\t%d\n',nY);
    fprintf(fid,'XStep\t%.6f\n',abs(dx));
    fprintf(fid,'YStep\t%.6f\n',abs(dy));
    fprintf(fid,'AcqE1\t0\nAcqE2\t0\nAcqE3\t0\n');
    fprintf(fid,['Euler angles refer to Sample Coordinate system (CS0)!\t' ...
                 'Mag\t100\tCoverage\t100\tDevice\t0\tKV\t20\t' ...
                 'TiltAngle\t70\tTiltAxis\t0\n']);
    fprintf(fid,'Phases\t%d\n',numel(phases));
    for k = 1:numel(phases)
        p = phases(k);
        fprintf(fid,'%.3f;%.3f;%.3f\t%.3f;%.3f;%.3f\t%s\t%d\t%d\n', ...
            p.abc(1),p.abc(2),p.abc(3), p.ang(1),p.ang(2),p.ang(3), ...
            p.name, p.laue, p.sg);
    end
    fprintf(fid,'Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n');

    M = [d.phase, d.x, d.y, d.bands, d.err, d.e1, d.e2, d.e3, d.mad, d.bc, d.bs].';
    fprintf(fid,'%d\t%.4f\t%.4f\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n',M);
end