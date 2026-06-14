function R = analyzeStrainRSM(rsm, options)
%ANALYZESTRAINRSM  Strain & composition of a film from an asymmetric RSM.
%   R = xrdc.rsm.analyzeStrainRSM(rsm, Substrate=, Film=, Reflection=, ...)
%
%   Auto-detects the substrate and film peaks in an asymmetric reciprocal-
%   space map, inverts them to in-plane (a_par) and out-of-plane (a_perp)
%   lattice parameters, decomposes the film's biaxial strain to a relaxed
%   pseudocubic a0, reports the degree of relaxation against the MEASURED
%   substrate in-plane parameter, and (for an alloy film with a Vegard table,
%   e.g. PZT) the composition x. Reports numbers for the DECLARED film; it
%   does not rank candidates.
%
%   1/d convention (xrdc.rsm.toReciprocalSpace): for reflection [h k l],
%   a_par = sqrt(h^2+k^2)/kPar, a_perp = l/kPerp.
%
%   rsm        : folder/file path (-> xrdc.rsm.loadAreaScan) or a pre-loaded
%                slice struct array (each as xrdc.io.emptyScan).
%   Substrate  : declared substrate material name (default "SrTiO3").
%   Film       : declared film material name (required; must be in materials.json).
%   Reflection : 1x3 [h k l]; an asymmetric reflection (l~=0 and h^2+k^2>0).
%   NoiseFactor : forwarded to xrdc.rsm.findRsmPeaks.
%
%   Output R: see field assignments below. R.elasticModel is "nu" or "c13c33";
%   R.nu is the Poisson ratio under the "nu" model and NaN under "c13c33"
%   (no single nu exists there -- the strain still uses the correct factor).
%   Caveats (pseudocubic a0; strain/composition real-data-unvalidated this
%   iteration) in SCIENTIFIC_ASSUMPTIONS.
%   Errors: xrdc:rsm:badReflection, xrdc:lattice:unknownMaterial.

    arguments
        rsm
        options.Substrate   (1,1) string  = "SrTiO3"
        options.Film        (1,1) string
        options.Reflection  (1,3) double
        options.NoiseFactor (1,1) double = 5
    end

    % Required name-value args (no default => field absent if omitted).
    if ~isfield(options, 'Reflection')
        error('xrdc:rsm:missingReflection', 'Reflection=[h k l] is required.');
    end
    if ~isfield(options, 'Film')
        error('xrdc:rsm:missingFilm', 'Film="<material>" is required.');
    end

    hkl = options.Reflection;
    hk2 = hkl(1)^2 + hkl(2)^2;
    if hkl(3) == 0 || hk2 == 0
        error('xrdc:rsm:badReflection', ...
            'Reflection [%g %g %g] is not asymmetric (need l~=0 and h^2+k^2>0).', ...
            hkl(1), hkl(2), hkl(3));
    end

    subMat  = xrdc.lattice.loadMaterials(options.Substrate);
    filmMat = xrdc.lattice.loadMaterials(options.Film);

    % --- load slices ---
    if isstruct(rsm)
        scans = rsm;
    elseif ischar(rsm) || isstring(rsm)
        rsmStr = string(rsm);
        if isfile(rsmStr)
            % Single file path: wrap in cell so loadAreaScan uses its file-list branch.
            scans = xrdc.rsm.loadAreaScan({char(rsmStr)});
        else
            scans = xrdc.rsm.loadAreaScan(rsmStr);
        end
    else
        scans = xrdc.rsm.loadAreaScan(rsm);
    end

    % --- assemble (kPar, kPerp, intensity) cloud ---
    kp = []; kz = []; I = [];
    for j = 1:numel(scans)
        [a, b] = xrdc.rsm.toReciprocalSpace(scans(j));
        kp = [kp; a(:)]; kz = [kz; b(:)]; I = [I; double(scans(j).counts(:))]; %#ok<AGROW>
    end

    pk = xrdc.rsm.findRsmPeaks(kp, kz, I, 'NoiseFactor', options.NoiseFactor);

    sqHK = sqrt(hk2);
    invPar  = @(kpar)  sqHK / abs(kpar);
    invPerp = @(kperp) hkl(3) / kperp;

    % --- substrate ---
    subAMeas = invPar(pk.substrate.kPar);
    subCMeas = invPerp(pk.substrate.kPerp);
    % prediction from declared lattice (uses tabulated a and c; for a cubic
    % substrate c=a, for a tetragonal substrate e.g. TiO2 they differ)
    kParPred  = sqHK / subMat.a;
    kPerpPred = hkl(3) / subMat.c;
    offPct = 100 * hypot(pk.substrate.kPar - kParPred, pk.substrate.kPerp - kPerpPred) ...
                 / hypot(kParPred, kPerpPred);
    flags = pk.flags;
    if offPct > 2
        flags(end+1) = "substrateOffPrediction"; %#ok<AGROW>
    end
    R.substrate = struct('name', string(subMat.name), 'found', pk.substrate.found, ...
        'kPar', pk.substrate.kPar, 'kPerp', pk.substrate.kPerp, ...
        'aMeas', subAMeas, 'cMeas', subCMeas, 'predictOffPct', offPct);

    % --- film ---
    if ~pk.film.found
        R.film = struct('found', false, 'kPar', NaN, 'kPerp', NaN);
        [R.aPar, R.aPerp, R.cFilm, R.a0, R.strainPar, R.strainPerp, ...
         R.relaxation, R.x] = deal(NaN);
        R.pseudomorphic = false;
        R.reflection = hkl; R.elasticModel = elasticModelName(filmMat);
        R.nu = nuOf(filmMat); R.flags = flags;
        R.notes = "film peak not found; substrate-only result";
        return
    end

    aPar  = invPar(pk.film.kPar);
    aPerp = invPerp(pk.film.kPerp);
    f = xrdc.lattice.elasticFactor(filmMat);
    [a0, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(aPar, aPerp, f);

    % relaxation vs MEASURED substrate in-plane parameter
    denom = a0 - subAMeas;
    if abs(denom) < 1e-6
        relax = NaN;
    else
        relax = (aPar - subAMeas) / denom;
    end
    pseudomorphic = isfinite(relax) && abs(relax) < 0.1;

    % composition (alloy films with a Vegard table only)
    x = NaN;
    if ~isempty(filmMat.composition)
        ax = filmMat.composition.a; xx = filmMat.composition.x;
        if a0 >= min(ax) && a0 <= max(ax)
            x = interp1(ax, xx, a0);
        else
            flags(end+1) = "compositionOutOfRange"; %#ok<AGROW>
        end
    end

    R.film = struct('found', true, 'kPar', pk.film.kPar, 'kPerp', pk.film.kPerp);
    R.aPar = aPar; R.aPerp = aPerp; R.cFilm = aPerp; R.a0 = a0;
    R.strainPar = epsPar; R.strainPerp = epsPerp;
    R.relaxation = relax; R.pseudomorphic = pseudomorphic; R.x = x;
    R.reflection = hkl; R.elasticModel = elasticModelName(filmMat);
    R.nu = nuOf(filmMat); R.flags = flags;
    R.notes = "a0 is a pseudocubic strain-model average (see SCIENTIFIC_ASSUMPTIONS)";
end

function name = elasticModelName(m)
    if isfield(m.elastic, 'nu') && ~isempty(m.elastic.nu)
        name = "nu";
    else
        name = "c13c33";
    end
end

function v = nuOf(m)
    if isfield(m.elastic, 'nu') && ~isempty(m.elastic.nu)
        v = m.elastic.nu;
    else
        v = NaN;     % c13/c33 model: no single nu
    end
end
