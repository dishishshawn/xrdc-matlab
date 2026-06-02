function doc = parseXmlFile(path)
%PARSEXMLFILE  Parse an XML file robustly for instrument exports.
%   doc = xrdc.io.parseXmlFile(path)
%
%   Uses JAXP DocumentBuilder configured to skip external-entity / DTD
%   resolution (avoids network fetches and CVE-2020-27201 class issues)
%   and to relax secure-processing limits that can trip on large
%   reciprocal-space-map XRDML files. Passing the file as java.io.File
%   bypasses xmlread's URI-resolution path, which hangs on paths
%   containing spaces or parentheses.
%
%   Input
%     path : path to an XML file (PANalytical XRDML, Bruker XML, etc.)
%
%   Output
%     doc  : org.w3c.dom.Document root, same shape returned by xmlread.

    arguments
        path (1,1) string
    end

    factory = javax.xml.parsers.DocumentBuilderFactory.newInstance();
    setFeatureIfSupported(factory, 'http://apache.org/xml/features/nonvalidating/load-external-dtd', false);
    setFeatureIfSupported(factory, 'http://xml.org/sax/features/external-general-entities', false);
    setFeatureIfSupported(factory, 'http://xml.org/sax/features/external-parameter-entities', false);

    builder = factory.newDocumentBuilder();

    % Silence the EntityResolver — anything <!ENTITY ...> referenced
    % outside the file gets an empty InputSource instead of a network fetch.
    builder.setEntityResolver(org.xml.sax.helpers.DefaultHandler());

    doc = builder.parse(java.io.File(char(path)));
end

function setFeatureIfSupported(factory, name, value)
    try
        factory.setFeature(name, value);
    catch
        % Feature not supported on this JAXP implementation — ignore.
    end
end
