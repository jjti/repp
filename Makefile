LOCAL_BIN=/usr/local/bin
APP=$(LOCAL_BIN)/defrag
APP_DATA=$${HOME}/.defrag
PRIMER3_CONF=$(APP_DATA)/primer3_config/

SETTINGS=./config/config.yaml
TEMP_SETTINGS=./config.temp.yaml

PLATFORM:=$(shell uname)

.PHONY: test

ifeq ($(PLATFORM),Windows_NT)
	$(error Windows not supported via make)
endif

all: build install

install:
	mkdir -p $(APP_DATA)

	# copy app and primer3 settings
	cp $(SETTINGS) $(TEMP_SETTINGS)
	echo primer3_config-path: $(PRIMER3_CONF) >> $(TEMP_SETTINGS)
	mv $(TEMP_SETTINGS) $(APP_DATA)/config.yaml
	cp -r ./vendor/primer3_config $(PRIMER3_CONF)

	# copy DNA, feature, and enzyme databases
	cp -r ./assets/addgene/db/** $(APP_DATA) 
	cp -r ./assets/igem/db/** $(APP_DATA)
	cp -r ./assets/dnasu/db/** $(APP_DATA)
	cp ./assets/snapgene/features.tsv $(APP_DATA)
	cp ./assets/neb/enzymes.tsv $(APP_DATA)

	chmod a+x ./bin/*

	# copy binaries to /usr/local/bin
ifeq ($(PLATFORM),Linux)
	cp ./bin/linux $(APP)
	chmod a+x ./vendor/linux/*
	cp -n ./vendor/linux/blastn $(LOCAL_BIN) || true
	cp -n ./vendor/linux/ntthal $(LOCAL_BIN) || true
	cp -n ./vendor/linux/primer3_core $(LOCAL_BIN) || true
	cp -n ./vendor/linux/blastdbcmd $(LOCAL_BIN) || true
endif

	# copy binaries to /usr/local/bin
ifeq ($(PLATFORM),Darwin)
	cp ./bin/darwin $(APP)
	chmod a+x ./vendor/darwin/*
	cp -n ./vendor/darwin/blastn $(LOCAL_BIN) || true
	cp -n ./vendor/darwin/ntthal $(LOCAL_BIN) || true
	cp -n ./vendor/darwin/primer3_core $(LOCAL_BIN) || true
	cp -n ./vendor/darwin/blastdbcmd $(LOCAL_BIN) || true
endif

build:
	go get
	env GOOS=linux go build -o ./bin/linux -v
	env GOOS=darwin go build -o ./bin/darwin -v
	env GOOS=windows go build -o ./bin/windows/defrag.exe -v

dbs:
	cd assets && sh makeblastdbs.sh

clean:
	rm -rf $(PRIMER3_CONF)
	rm $(APP_DATA)/*

uninstall: clean
	rm $(APP)
	rm -rf $(APP_DATA)

test:
	go test ./internal/defrag

e2e:
	go test ./internal/defrag -tags=e2e