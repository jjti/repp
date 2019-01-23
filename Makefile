# Basic go commands
GOCMD=go
GOBUILD=$(GOCMD) build
GOCLEAN=$(GOCMD) clean
GOTEST=$(GOCMD) test
GOGET=$(GOCMD) get

DEFRAGBIN=/usr/local/bin/defrag
DEFRAGCONF=$${HOME}/.defrag

PLATFORM:=$(shell uname)

# copy settings file, going to add to it during install
SETTINGS=./config/config.yaml
TEMPSETTINGS=./defrag.temp.yaml

ifeq ($(OS),Windows_NT)
	$(error "defrag" does not support Windows)
endif

install: build
	mkdir -p $(DEFRAGCONF)
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
	brew install blast # install BLAST
endif
ifeq (, $(shell which primer3_core))
	brew install primer3 # install primer3
endif
	echo "\nprimer3_config-path: /usr/local/share/primer3/primer3_config/" >> $(TEMPSETTINGS)
	cp ./bin/darwin $(DEFRAGBIN)
endif

	# copy config file to /etc/defrag
	mv $(TEMPSETTINGS) $(DEFRAGCONF)/config.yaml

	# copy BLAST databases
	cp -r ./assets/addgene/db/** $(DEFRAGCONF) 
	cp -r ./assets/igem/db/** $(DEFRAGCONF)
	
build:
	# build for supported operating systems
	env GOOS=linux $(GOBUILD) -o ./bin/linux -v
	env GOOS=darwin $(GOBUILD) -o ./bin/darwin -v

uninstall:
	# removing defrag from filesystem
	rm $(DEFRAGBIN)
	rm -rf $(DEFRAGCONF)
