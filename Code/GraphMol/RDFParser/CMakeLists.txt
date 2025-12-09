rdkit_library(RDFParser
        RDFParser.cpp
        LINK_LIBRARIES GenericGroups Depictor GraphMol ChemReactions )
target_compile_definitions(RDFParser PRIVATE RDKIT_RDFPARSER_BUILD)

rdkit_headers(RDFParser.h
        DEST GraphMol/RDFParser)

rdkit_catch_test(testRDFtoMol testRDFtoMol.cpp
        LINK_LIBRARIES RDFParser FileParsers SubstructMatch )
