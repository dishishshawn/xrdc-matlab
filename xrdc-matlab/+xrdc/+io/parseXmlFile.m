function doc = parseXmlFile(path)
%PARSEXMLFILE  Parse an XML file robustly for instrument exports.
%   doc = xrdc.io.parseXmlFile(path)
%
%   Uses matlab.io.xml.dom (MATLAB's native, Java-free DOM). This is the
%   deployment-critical choice: the legacy javax.xml / JAXP path is NOT
%   available in a compiled standalone (the MATLAB Runtime has no JVM on
%   the resolution path), so the old DocumentBuilderFactory call failed in
%   XRDC.exe with "Unable to resolve the name
%   'javax.xml.parsers.DocumentBuilderFactory.newInstance'". The native
%   parser also takes a plain filename, sidestepping xmlread's URI
%   resolution that hangs on paths with spaces or parentheses, and does not
%   resolve external entities/DTDs by default (no network fetch, no XXE) —
%   the same hardening the JAXP version configured by hand.
%
%   Input
%     path : path to an XML file (PANalytical XRDML, Bruker XML, etc.)
%
%   Output
%     doc  : matlab.io.xml.dom.Document. W3C-DOM-compatible API used by the
%            readers: getDocumentElement, getElementsByTagName, item,
%            getLength, hasAttribute, getAttribute, getTextContent.

    arguments
        path (1,1) string
    end

    % The native parser is secure by default: it does not resolve external
    % entities/DTDs over the network, so no extra hardening is needed (unlike
    % the old JAXP path, which had to disable those features by hand).
    parser = matlab.io.xml.dom.Parser;
    doc = parseFile(parser, char(path));
end
