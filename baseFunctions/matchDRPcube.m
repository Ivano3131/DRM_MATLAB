function res = matchDRPcube(drpIn, drpDic, euDic, exp_para, options)
% matchDRPcube  Match measured DRPs to the cube-lath dictionary.
%
%   res = matchDRPcube(drp_original, drpDic, euDic, exp_para)
%   res = matchDRPcube(oneDRP,       drpDic, euDic, exp_para)
%
% Brute-force normalised cross-correlation over the whole dictionary AND all
% ph_num azimuth shifts.  No autoencoder, no kNN, no latent space - the score
% is the actual correlation between the measured and the simulated pattern,
% so it is interpretable (1 = identical, 0 = unrelated) and there is no
% training step to go wrong.
%
% HOW phi1 COMES OUT ------------------------------------------------------
% The dictionary is built at phi1 = 0 (makeDRPdic_cube).  Increasing phi1
% rotates the whole DRP in azimuth, which on this grid is exactly
% circshift(sim, s, 2) with phi1 = s * 360/ph_num.  So
%
%       score(k,s) = < meas , circshift(sim_k, s, 2) >
%
% after zero-meaning and unit-normalising both patterns, and the best (k,s)
% gives Phi and phi2 from euDic(k,:) and phi1 from s.  All ph_num shifts are
% evaluated at once by the circular correlation theorem,
%
%       score(k,:) = real( ifft( sum_theta fft(meas) .* conj(fft(sim_k)) ) )
%
% so the shift search costs one FFT, not ph_num inner products.  This also
% removes the failure mode of the older path, where the dictionary was
% canonicalised in azimuth by column sum (makeDRPdic_Ti64.m:61) but the
% measurement by count-of-max-pixels-per-column (IndexingEngine.m:37) - two
% different rules, biasing phi1.
%
% INPUTS ------------------------------------------------------------------
%   drpIn     n1 x n2 cell of th_num x ph_num DRPs (from igray2drp), or a
%             single th_num x ph_num numeric DRP.
%   drpDic    n x 1 cell from makeDRPdic_cube.
%   euDic     n x 2 [Phi phi2] in degrees.
%   exp_para  needs th_num, ph_num.
%
% OPTIONS -----------------------------------------------------------------
%   chunk    pixels processed per batch.  0 (default) picks it automatically:
%            32 on the CPU, and on the GPU as much as free device memory
%            allows.  Pixels are independent columns of the batched product,
%            so chunk changes speed and memory only, never the result.
%   gpu      "auto" (default) uses the GPU when Parallel Computing Toolbox and
%            a supported device are present, and falls back to the CPU
%            silently otherwise.  "on" errors if there is no device, "off"
%            forces the CPU path.
%   verbose  true (default) prints progress.
%
% GPU NOTES ---------------------------------------------------------------
% The device path runs in DOUBLE, like the CPU path, so the two agree to
% round-off (~1e-15 on a score in [-1,1]) rather than to single precision.
% That matters because the winner is picked by argmax over nDic*ph_num
% candidates: a coarser precision could reorder near-ties.  Consumer cards
% run FP64 at a fraction of FP32, so the win here comes mostly from
% bandwidth - the ifft/max stage moves ~2.6 MB per pixel and is what the
% device actually accelerates.
%
% Note the dictionary is degenerate at Phi = 0, where every phi2 gives the
% same pattern (phi1 absorbs the rotation about c).  Those entries tie
% exactly, so CPU and GPU may report different .idx / phi2 for them while
% describing the identical DRP and the identical orientation.  Compare
% orientations or .score there, not .idx.
%
% OUTPUT ------------------------------------------------------------------
%   res.EUmap  n1 x n2 x 3 [phi1 Phi phi2] in degrees
%   res.score  n1 x n2 correlation in [-1,1].  HIGHER IS BETTER - note this
%              is the opposite sense to IndexingEngine's .quality, which is a
%              distance.
%   res.shift  n1 x n2 azimuth shift in ph_num steps
%   res.idx    n1 x n2 dictionary index
%   res.euler  n1*n2 x 3, the same orientations flattened (handy for a single
%              DRP, where it is just 1 x 3)
%
% See also makeDRPdic_cube, DRPsim_cube, check_indexing_result_cube.
% -------------------------------------------------------------------------
arguments
    drpIn
    drpDic cell
    euDic double
    exp_para struct
    options.chunk (1,1) double {mustBeNonnegative} = 0
    options.gpu (1,1) string {mustBeMember(options.gpu,["auto","on","off"])} = "auto"
    options.verbose (1,1) logical = true
end

th_num = exp_para.th_num;
ph_num = exp_para.ph_num;

% ---- accept a single DRP as well as a full map --------------------------
if ~iscell(drpIn)
    drpCell = {drpIn};
else
    drpCell = drpIn;
end
[n1,n2] = size(drpCell);
P       = n1 * n2;
nDic    = numel(drpDic);

if size(euDic,1) ~= nDic
    error('matchDRPcube:euDic','euDic must have one row per dictionary entry.');
end

% ---- pre-normalise and pre-transform the dictionary --------------------
% Zero mean + unit Frobenius norm, so the correlation is bounded by 1 and is
% invariant to the brightness and offset of the measurement.  A circular
% shift changes neither, so this survives the shift search.
D = zeros(th_num, ph_num, nDic);
for k = 1:nDic
    D(:,:,k) = local_normalise(drpDic{k});
end
% conj(FFT) along azimuth, laid out as nDic x th_num x ph_num so that each
% azimuth FREQUENCY is one page of a batched matrix product
A = permute(conj(fft(D,[],2)), [3 1 2]);
clear D

useFast = exist('pagemtimes','builtin') || exist('pagemtimes','file');

% ---- device selection ---------------------------------------------------
[useGPU, gpuName, gpuFree] = local_pickDevice(options.gpu);

% Below a batch or so of pixels the dictionary transfer dominates and the
% device is a net loss - cubeAmbiguity calls this once per test pattern.
if useGPU && P < 64 && options.gpu ~= "on"
    useGPU = false;
end

% ---- batch size --------------------------------------------------------
% Each pixel needs the complex product (nDic x ph_num x 16 B), its real
% inverse transform (x 8 B) and headroom for the transform's own workspace.
chunkSz = options.chunk;
if chunkSz == 0
    if useGPU
        bytesPerPix = nDic * ph_num * 40;
        chunkSz = max(64, min(2048, floor(0.5 * gpuFree / bytesPerPix)));
    else
        chunkSz = 32;
    end
end
chunkSz = min(chunkSz, P);

Acpu = [];
if useGPU
    Acpu = A;                      % kept for the out-of-memory retry below
    A    = gpuArray(A);
end

score = zeros(P,1);
shift = zeros(P,1);
idx   = ones(P,1);

if options.verbose
    fprintf('matchDRPcube: %d pixels against %d patterns x %d shifts ...\n', ...
        P, nDic, ph_num);
    if useGPU
        fprintf('  device: %s (double), chunk %d\n', gpuName, chunkSz);
    else
        fprintf('  device: CPU (double), chunk %d\n', chunkSz);
    end
end
tStart = tic;

i0     = 1;
nBatch = 0;
while i0 <= P
    i1 = min(i0 + chunkSz - 1, P);
    pIdx = i0:i1;
    nP   = numel(pIdx);

    % measured block, normalised the same way
    Mblk  = zeros(th_num, ph_num, nP);
    valid = false(nP,1);
    for jj = 1:nP
        [Mblk(:,:,jj), valid(jj)] = local_normalise(drpCell{pIdx(jj)});
    end

    if useGPU
        try
            [best, kBest, sBest] = local_matchGPU(A, Mblk, nP, ph_num);
        catch ME
            % A device that cannot hold this batch is a reason to use a
            % smaller one, not to lose the run.
            if contains(ME.identifier,'OutOfMemory') && chunkSz > 64
                chunkSz = max(64, floor(chunkSz/2));
                if options.verbose
                    fprintf('  GPU out of memory, retrying with chunk %d\n', chunkSz);
                end
                reset(gpuDevice);          % drops every gpuArray, A included
                A = gpuArray(Acpu);
                continue                       % same i0, smaller batch
            end
            rethrow(ME);
        end
    else
        % th_num x nP x ph_num, pages over azimuth frequency
        Fm = permute(fft(Mblk,[],2), [1 3 2]);

        if useFast
            C = pagemtimes(A, Fm);                 % nDic x nP x ph_num
        else
            C = complex(zeros(nDic, nP, ph_num));
            for f = 1:ph_num
                C(:,:,f) = A(:,:,f) * Fm(:,:,f);
            end
        end

        % circular correlation theorem: index 1 along dim 3 is shift 0
        cc = real(ifft(C, [], 3));
        cc = reshape(permute(cc, [1 3 2]), nDic*ph_num, nP);

        [best, li] = max(cc, [], 1);
        kBest = mod(li-1, nDic) + 1;
        sBest = floor((li-1) / nDic);
    end

    best(~valid)  = 0;
    kBest(~valid) = 1;
    sBest(~valid) = 0;

    score(pIdx) = best(:);
    idx(pIdx)   = kBest(:);
    shift(pIdx) = sBest(:);

    if options.verbose && (i1 == P || mod(nBatch, 50) == 0)
        fprintf('  %6d / %6d pixels  (%.1f s)\n', i1, P, toc(tStart));
    end
    nBatch = nBatch + 1;
    i0 = i1 + 1;
end

if useGPU
    A = [];                                    %#ok<NASGU> release device memory
end

if options.verbose
    fprintf('matchDRPcube: done in %.1f s, median score %.3f.\n', ...
        toc(tStart), median(score));
end

% ---- assemble ----------------------------------------------------------
phi1 = mod(shift * (360/ph_num), 360);
euler = [phi1, euDic(idx,1), euDic(idx,2)];
euler(score == 0, :) = 0;      % pixels with no signal (flat/blank DRPs)

res.EUmap = reshape(euler, n1, n2, 3);
res.score = reshape(score, n1, n2);
res.shift = reshape(shift, n1, n2);
res.idx   = reshape(idx,   n1, n2);
res.euler = euler;
end

% =========================================================================
function [best, kBest, sBest] = local_matchGPU(A, Mblk, nP, ph_num)
% One batch on the device.  Same arithmetic as the CPU branch, in double, so
% the two agree to round-off.
Fm = permute(fft(gpuArray(Mblk), [], 2), [1 3 2]);   % th_num x nP x ph_num
C  = pagemtimes(A, Fm);                              % nDic   x nP x ph_num
clear Fm
cc = real(ifft(C, [], 3));
clear C

% Argmax in two stages rather than one over a reshaped nDic*ph_num x nP
% copy.  The CPU branch's linear index is k + s*nDic, so its max takes the
% smallest s, and within it the smallest k.  max over dim 1 keeps the
% smallest k per shift and max over dim 3 then keeps the smallest shift -
% the same winner, without materialising the permute (which on this array is
% another 0.66 MB per pixel of device traffic).
[m1, kAll]   = max(cc, [], 1);          % 1 x nP x ph_num
[best, sIdx] = max(m1, [], 3);          % 1 x nP
kBest = reshape(kAll, nP, ph_num);
kBest = kBest((1:nP).' + (sIdx(:)-1)*nP).';
sBest = sIdx - 1;

[best, kBest, sBest] = gather(best, kBest, sBest);
end

% =========================================================================
function [useGPU, name, freeBytes] = local_pickDevice(mode)
% "auto" never costs the caller a run: anything missing means the CPU path.
useGPU    = false;
name      = '';
freeBytes = 0;
if mode == "off"
    return
end
why = '';
if isempty(ver('parallel'))
    why = 'Parallel Computing Toolbox is not installed';
elseif gpuDeviceCount("available") < 1
    why = 'no supported GPU device is available';
else
    try
        d = gpuDevice;
        if ~d.DeviceAvailable
            why = 'the GPU device is not available for computation';
        else
            useGPU    = true;
            name      = d.Name;
            freeBytes = d.AvailableMemory;
        end
    catch ME
        why = ME.message;
    end
end
if ~useGPU && mode == "on"
    error('matchDRPcube:noGPU','gpu="on" was requested but %s.', why);
end
end

% =========================================================================
function [d, ok] = local_normalise(d)
% Zero mean, unit Frobenius norm.  ok = false for a flat pattern (background),
% which carries no orientation information.
d = double(d);
d = d - mean(d(:));
nrm = norm(d, 'fro');
ok = nrm > 0;
if ok
    d = d / nrm;
else
    d = zeros(size(d));
end
end
