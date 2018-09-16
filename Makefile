# Basic go commands
GOCMD=/usr/local/go/bin/go
GOBUILD=$(GOCMD) build
GOCLEAN=$(GOCMD) clean
GOTEST=$(GOCMD) test
GOGET=$(GOCMD) get

# Binary name
BINARY_NAME=./bin/decvec

# Tool directory
TOOL_DIR=./tools/

# BLAST paths
# TODO: for different operating systems
BLAST_NAME=ncbi-blast-2.7.1+
BLAST_MAC_FTP=ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/$(BLAST_NAME)-x64-macosx.tar.gz
BLAST_MAC_LOCAL=$(TOOL_DIR)$(BLAST_NAME).tar.gz

# Primer3 paths
PRIMER3_NAME=primer3-2.4.0
PRIMER3_SOURCEFORGE=https://sourceforge.net/projects/primer3/files/primer3/2.4.0/primer3-2.4.0.tar.gz/download
PRIMER3_LOCAL=$(TOOL_DIR)$(PRIMER3_NAME).tar.gz

all: build mock
build: 
		$(GOBUILD) -o $(BINARY_NAME) -v
test: 
		$(GOTEST) -v ./...
clean: 
		$(GOCLEAN)
		rm -f $(BINARY_NAME)
		rm -f $(BINARY_UNIX)
run:
		$(GOBUILD) -o $(BINARY_NAME) -v ./...
		./$(BINARY_NAME)
deps: 
		# download BLAST
		curl $(BLAST_MAC_FTP) -o $(BLAST_MAC_LOCAL)
		tar xvzf $(BLAST_MAC_LOCAL) -C $(TOOL_DIR)
		rm $(BLAST_MAC_LOCAL)

		# download primer3
		curl -L $(PRIMER3_SOURCEFORGE) -o $(PRIMER3_LOCAL)
		tar xvzf $(PRIMER3_LOCAL) -C $(TOOL_DIR)
		rm $(PRIMER3_LOCAL)

		# make primer3
		make -C $(TOOL_DIR)$(PRIMER3_NAME)/src
mock: 
		$(BINARY_NAME) make -t ../test/payloads/target.fa -f ../test/payloads/input.fa
