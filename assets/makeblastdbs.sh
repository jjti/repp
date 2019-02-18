#!/bin/sh

cd addgene/db
makeblastdb -in ./addgene -parse_seqids -title addgene -dbtype nucl

cd ../../igem/db
makeblastdb -in ./igem -parse_seqids -title igem -dbtype nucl

cd ../../dnasu/db
makeblastdb -in ./dnasu -parse_seqids -title dnasu -dbtype nucl