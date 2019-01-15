# Basic go commands
GOCMD=go
GOBUILD=$(GOCMD) build
GOCLEAN=$(GOCMD) clean
GOTEST=$(GOCMD) test
GOGET=$(GOCMD) get

DEFRAGBIN=/usr/local/bin/defrag
DEFRAGETC=/etc/defrag

PLATFORM := $(shell uname)

# copy settings file, going to add to it during install
SETTINGS=./settings.yaml
TEMPSETTINGS=./settings.temp

ifeq ($(OS),Windows_NT)
	$(error "defrag" does not support Windows)
endif

install:
	mkdir -p $(DEFRAGETC)
	cp $(SETTINGS) $(TEMPSETTINGS)

# install outside dependencies, copy binary to /usr/local/bin
ifeq ($(PLATFORM),Linux)
	apt-get install ncbi-blast+ primer3

	echo "\nprimer3_config-path: /etc/primer3_config/" >> $(TEMPSETTINGS)

	cp ./bin/linux $(DEFRAGBIN)
endif

ifeq ($(PLATFORM),Darwin)
ifeq (, $(shell which brew))
	# install homebrew
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
endif
ifeq (, $(shell which blastn))
	# install BLAST
	brew install blast
endif
ifeq (, $(shell which primer3_core))
	# install primer3
	brew install primer3
endif
	echo "\nprimer3_config-path: /usr/local/share/primer3/primer3_config/" >> $(TEMPSETTINGS)

	cp ./bin/darwin $(DEFRAGBIN)
endif

	# install config file to /etc/defrag
	mv $(TEMPSETTINGS) $(DEFRAGETC)/settings.yaml

	# unzip BLAST databases
ifeq (,$(wildcard ./addgene.zip))
		unzip -jo ./addgene.zip -d $(DEFRAGETC) 
		unzip -jo ./igem.zip -d $(DEFRAGETC) 
endif
	
build:
	# build for all operating systems
	env GOOS=linux $(GOBUILD) -o ./bin/linux -v
	env GOOS=darwin $(GOBUILD) -o ./bin/darwin -v

uninstall:
	# removing defrag from filesystem
	rm $(DEFRAGBIN)
	rm -rf $(DEFRAGETC)


