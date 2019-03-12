GOCMD=go

APP=/usr/local/bin/defrag
APPDATA=$${HOME}/.defrag

SETTINGS=./config/config.yaml
TEMPSETTINGS=./config.temp.yaml
PLATFORM:=$(shell uname)

.PHONY: test

ifeq ($(PLATFORM),Windows_NT)
	$(error "defrag" does not support Windows)
endif

install:
	mkdir -p $(APPDATA)
	cp $(SETTINGS) $(TEMPSETTINGS)

ifeq ($(PLATFORM),Linux)
	apt-get update
	apt-get install ncbi-blast+ primer3
	echo "\nprimer3_config-path: /etc/primer3_config/" >> $(TEMPSETTINGS)
	cp ./bin/linux $(APP)
endif

ifeq ($(PLATFORM),Darwin)
ifeq (, $(shell which brew))
	/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
endif
ifeq (, $(shell which blastn))
	brew install blast
endif
ifeq (, $(shell which primer3_core))
	brew install primer3
endif
	echo "\nprimer3_config-path: /usr/local/share/primer3/primer3_config/" >> $(TEMPSETTINGS)
	cp ./bin/darwin $(APP)
endif

	mv $(TEMPSETTINGS) $(APPDATA)/config.yaml

	cp -r ./assets/addgene/db/** $(APPDATA) 
	cp -r ./assets/igem/db/** $(APPDATA)
	cp -r ./assets/dnasu/db/** $(APPDATA)
	cp ./assets/snapgene/features.tsv $(APPDATA)
	cp ./assets/neb/enzymes.tsv $(APPDATA)
	
build:
	go get
	env GOOS=linux go build -o ./bin/linux -v
	env GOOS=darwin go build -o ./bin/darwin -v

dbs:
	cd assets && sh makeblastdbs.sh

clean:
	rm $(APPDATA)/*

uninstall: clean
	rm $(APP)
	rm -rf $(APPDATA)

test:
	go test ./internal/defrag

e2e:
	go test ./internal/defrag -tags=e2e