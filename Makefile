# Basic go commands
GOCMD=go
GOBUILD=$(GOCMD) build
GOCLEAN=$(GOCMD) clean
GOTEST=$(GOCMD) test
GOGET=$(GOCMD) get

ETC_DIR=/etc/defrag

DEFRAG_DBS=https://console.cloud.google.com/storage/browser/defrag

# get the platform type
PLATFORM := $(shell uname)

# copy settings file, going to add to it during install
SETTINGS=./settings.yaml
TEMP_SETTINGS=./settings.temp

ifeq ($(OS),Windows_NT)
	$(error "defrag" does not support Windows)
endif

install:
	mkdir -p $(ETC_DIR)
	cp $(SETTINGS) $(TEMP_SETTINGS)

# install outside dependencies, copy binary to /usr/local/bin
ifeq ($(PLATFORM),Linux)
	apt-get install ncbi-blast+ primer3

	echo "\nprimer3_config-path: /etc/primer3_config/" >> $(TEMP_SETTINGS)

	cp ./bin/linux /usr/local/bin/defrag
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
	echo "\nprimer3_config-path: /usr/local/share/primer3/primer3_config/" >> $(TEMP_SETTINGS)

	cp ./bin/darwin /usr/local/bin/defrag
endif

	# install config file to /etc/defrag
	mv $(TEMP_SETTINGS) $(ETC_DIR)/settings.yaml

	# unzip BLAST databases
	unzip -jo ./assets/addgene/addgene.zip -d $(ETC_DIR) 
	unzip -jo ./assets/igem/igem.zip -d $(ETC_DIR) 
	
build:
		# build for all operating systems
		env GOOS=linux $(GOBUILD) -o ./bin/linux -v
		env GOOS=darwin $(GOBUILD) -o ./bin/darwin -v


