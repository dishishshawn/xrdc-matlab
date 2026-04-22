function tests = testLattice
%TESTLATTICE  Unit tests for the +xrdc.+lattice package.
    tests = functiontests(localfunctions);
end

% ---------- Bragg's law round-trip ----------

function testBraggRoundTrip(tc)
    lambda = 1.5406;
    twoTheta = [20, 45, 90, 120];
    d = xrdc.lattice.twoThetaToD(twoTheta, lambda);
    twoThetaBack = xrdc.lattice.dToTwoTheta(d, lambda);
    tc.verifyEqual(twoThetaBack, twoTheta, 'AbsTol', 1e-10);
end

function testBraggKnownPeak(tc)
    % Si (400) with Cu Kalpha1 is tabulated at 2θ = 69.130°
    % d = a / √16 for cubic Si, a = 5.43088 Å → d = 1.35772 Å
    d = 5.43088 / 4;
    twoTheta = xrdc.lattice.dToTwoTheta(d, 1.5406);
    tc.verifyEqual(twoTheta, 69.13, 'AbsTol', 0.05);
end

function testBraggForbiddenReturnsNaN(tc)
    % λ = 1.54 Å, d = 0.5 Å ⇒ λ/(2d) = 1.54 > 1
    twoTheta = xrdc.lattice.dToTwoTheta(0.5, 1.54);
    tc.verifyTrue(isnan(twoTheta));
end

% ---------- Energy ↔ wavelength ----------

function testEnergyWavelengthRoundTrip(tc)
    % Cu Kα1 default in XRAY.def: 8049.19 eV
    lambda = xrdc.lattice.energyToLambda(8049.19);
    tc.verifyEqual(lambda, 1.5406, 'AbsTol', 5e-4);
    energyBack = xrdc.lattice.lambdaToEnergy(lambda);
    tc.verifyEqual(energyBack, 8049.19, 'RelTol', 1e-10);
end

% ---------- d-spacing for each crystal system ----------

function testCubicSi(tc)
    lat = struct('system','cubic','a', 5.43088);
    d = xrdc.lattice.dSpacingFromHKL(4, 0, 0, lat);
    tc.verifyEqual(d, 5.43088/4, 'AbsTol', 1e-10);
end

function testTetragonal(tc)
    % Known rutile (TiO2): a=4.593, c=2.959, (110) at d=3.2474 Å
    lat = struct('system','tetragonal','a',4.593,'c',2.959);
    d = xrdc.lattice.dSpacingFromHKL(1,1,0,lat);
    tc.verifyEqual(d, 4.593/sqrt(2), 'AbsTol', 1e-6);
end

function testOrthorhombic(tc)
    lat = struct('system','orthorhombic','a',5,'b',6,'c',7);
    % 1/d² = 1/25 + 1/36 + 1/49 for (1,1,1)
    dExpected = 1/sqrt(1/25 + 1/36 + 1/49);
    d = xrdc.lattice.dSpacingFromHKL(1,1,1,lat);
    tc.verifyEqual(d, dExpected, 'AbsTol', 1e-10);
end

function testHexagonal(tc)
    % GaN: a=3.189, c=5.185; (00.2) → d = c/2 = 2.5925
    lat = struct('system','hexagonal','a',3.189,'c',5.185);
    d = xrdc.lattice.dSpacingFromHKL(0,0,2,lat);
    tc.verifyEqual(d, 5.185/2, 'AbsTol', 1e-10);
end

function testMonoclinic(tc)
    % When β=90° the monoclinic formula should reduce to orthorhombic
    lat = struct('system','monoclinic','a',5,'b',6,'c',7,'beta',90);
    d = xrdc.lattice.dSpacingFromHKL(1,0,0,lat);
    tc.verifyEqual(d, 5, 'AbsTol', 1e-10);
    d = xrdc.lattice.dSpacingFromHKL(0,0,1,lat);
    tc.verifyEqual(d, 7, 'AbsTol', 1e-10);
end

function testTriclinicReducesToCubic(tc)
    % Triclinic with a=b=c, α=β=γ=90 should match cubic
    lat = struct('system','triclinic','a',4,'b',4,'c',4,...
                 'alpha',90,'beta',90,'gamma',90);
    d = xrdc.lattice.dSpacingFromHKL(2,0,0,lat);
    tc.verifyEqual(d, 2, 'AbsTol', 1e-10);
end

function testRhombohedralReducesToCubic(tc)
    % Rhombohedral with α=90° should match cubic
    lat = struct('system','rhombohedral','a',4,'alpha',90);
    d = xrdc.lattice.dSpacingFromHKL(1,0,0,lat);
    tc.verifyEqual(d, 4, 'AbsTol', 1e-10);
end

function testHKLVectorised(tc)
    lat = struct('system','cubic','a', 4);
    d = xrdc.lattice.dSpacingFromHKL([1;2;3], [0;0;0], [0;0;0], lat);
    tc.verifyEqual(d, [4; 2; 4/3], 'AbsTol', 1e-10);
end

% ---------- Nelson–Riley ----------

function testNelsonRileyRecovery(tc)
    % Construct a perfect NR line: a_i = 5 + 0.1 * NR(θ_i)
    twoTheta = [30, 50, 70, 90, 110, 130];
    theta    = deg2rad(twoTheta / 2);
    nr = cos(theta).^2 ./ sin(theta) + cos(theta).^2 ./ (twoTheta/2);
    a = 5 + 0.1 * nr;
    res = xrdc.lattice.nelsonRiley(twoTheta.', a.');
    tc.verifyEqual(res.a0, 5, 'AbsTol', 1e-10);
    tc.verifyEqual(res.slope, 0.1, 'AbsTol', 1e-10);
    tc.verifyEqual(res.rSquared, 1, 'AbsTol', 1e-10);
end

function testNelsonRileyWithNoise(tc)
    % Perturbation should not shift the intercept by much
    rng(42);
    twoTheta = linspace(50, 130, 6).';
    theta    = deg2rad(twoTheta / 2);
    nr = cos(theta).^2 ./ sin(theta) + cos(theta).^2 ./ (twoTheta/2);
    a = 5.43 + 0.01 * nr + 0.0005 * randn(size(nr));
    res = xrdc.lattice.nelsonRiley(twoTheta, a);
    tc.verifyEqual(res.a0, 5.43, 'AbsTol', 0.002);
    tc.verifyTrue(res.rSquared > 0.8);
    tc.verifyTrue(isfinite(res.a0SE) && res.a0SE > 0);
end

function testNelsonRileySEvsFitlm(tc)
    % Cross-check a0SE / slopeSE against fitlm (textbook OLS SEs).
    % Gated on the Statistics & Machine Learning Toolbox being present —
    % this guards against regressions that would reintroduce the Delphi
    % 1/√n scaling bug in xrdc3.pas:281.
    if isempty(which('fitlm'))
        tc.assumeFail('fitlm not available (Statistics Toolbox missing).');
    end
    rng(7);
    twoTheta = linspace(40, 130, 8).';
    theta    = deg2rad(twoTheta / 2);
    nr = cos(theta).^2 ./ sin(theta) + cos(theta).^2 ./ (twoTheta/2);
    a  = 5.43 + 0.008 * nr + 0.001 * randn(size(nr));

    res = xrdc.lattice.nelsonRiley(twoTheta, a);
    lm  = fitlm(nr, a);                           % y = b0 + b1 * x
    seIntercept = lm.Coefficients.SE(1);
    seSlope     = lm.Coefficients.SE(2);

    tc.verifyEqual(res.a0SE,    seIntercept, 'RelTol', 1e-8);
    tc.verifyEqual(res.slopeSE, seSlope,     'RelTol', 1e-8);
end

% ---------- Kiessig thickness ----------

function testKiessigRoundTrip(tc)
    % Generate synthetic fringes for a known thickness, then recover it
    lambda  = 1.5406;                     % Å
    t_nm    = 50;                          % true thickness 50 nm
    lambda_nm = lambda / 10;
    % sin θ_i = sin θ_0 + i * λ/(2t)
    slope = lambda_nm / (2 * t_nm);
    sinTh = 0.02 + (0:9).' * slope;
    twoTheta = 2 * rad2deg(asin(sinTh));
    res = xrdc.lattice.thicknessFromFringes(twoTheta, lambda);
    tc.verifyEqual(res.thicknessNm,     t_nm, 'AbsTol', 0.01);
    tc.verifyEqual(res.thicknessFitNm,  t_nm, 'AbsTol', 0.01);
end

% ---------- simulatePattern ----------

function testSimulateSrTiO3(tc)
    % SrTiO3 is cubic, a = 3.905 Å. (100), (200), (300), (400) at
    % 22.8, 46.48, 72.5, 104.32 deg (from Substrates.def / tabulated Cu Kα).
    lat = struct('system','cubic','a',3.905);
    T = xrdc.lattice.simulatePattern(lat, 4, 1.5406, ...
        'TwoThetaRange', [20, 110]);
    % extract the (h,0,0) families (symmetry-collapsed: the survivor is the
    % one with min h+k+l then min l*10000+k*100+h — which for {100} is (0,0,1))
    % So we just pick the rows with h²+k²+l² = 1, 4, 9, 16.
    sumSq = T.h.^2 + T.k.^2 + T.l.^2;
    d100 = T.twoTheta(sumSq == 1);  d100 = d100(1);
    d200 = T.twoTheta(sumSq == 4);  d200 = d200(1);
    d300 = T.twoTheta(sumSq == 9);  d300 = d300(1);
    d400 = T.twoTheta(sumSq == 16); d400 = d400(1);
    tc.verifyEqual(d100, 22.8,   'AbsTol', 0.05);
    tc.verifyEqual(d200, 46.48,  'AbsTol', 0.05);
    tc.verifyEqual(d300, 72.5,   'AbsTol', 0.1);
    tc.verifyEqual(d400, 104.32, 'AbsTol', 0.2);
end

function testSimulateDuplicatesCollapsed(tc)
    % In cubic, (1,0,0), (0,1,0), (0,0,1) all have the same d; they
    % should collapse to exactly one row when KeepEquivalents=false.
    lat = struct('system','cubic','a',4);
    T = xrdc.lattice.simulatePattern(lat, 1, 1.5406, ...
        'TwoThetaRange', [0, 180]);
    % Count rows where d ≈ 4
    nD4 = sum(abs(T.d - 4) < 1e-6);
    tc.verifyEqual(nD4, 1);
end
