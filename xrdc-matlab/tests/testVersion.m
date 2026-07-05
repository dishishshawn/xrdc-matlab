function tests = testVersion
%TESTVERSION  Unit tests for xrdc.version.
    tests = functiontests(localfunctions);
end

function testKnownAnswer(tc)
    tc.verifyEqual(xrdc.version(), "1.1.0");
end

function testSemverFormat(tc)
    v = xrdc.version();
    tc.verifyTrue(isstring(v) && isscalar(v));
    tc.verifyMatches(v, "^\d+\.\d+\.\d+$");
end
