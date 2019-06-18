LOCAL_BIN=/usr/local/bin
APP=${LOCAL_BIN}/plade
APP_DATA=$${HOME}/.plade
SETTINGS=./config/config.yaml

NAME=plasmid-defragger
VERSION=0.1.0

DIST=./dist
DIST_WIN_ZIP=${DIST}/${NAME}_windows_${VERSION}.zip
DIST_SRC=${DIST}/${NAME}_src_${VERSION}
DIST_SRC_TAR=${DIST_SRC}.tar.gz

PLATFORM:=$(shell uname)

.PHONY: test dist
.DEFAULT_GOAL: build

ifeq ($(PLATFORM),Windows_NT)
	$(error Windows not supported via make)
endif

build:
	rm -f ./bin/* 
	go get -d
	env GOOS=linux go build -o ./bin/linux -v
	env GOOS=darwin go build -o ./bin/darwin -v
	env GOOS=windows go build -o ./bin/plade.exe -v

install:
	mkdir -p $(APP_DATA)

	cp $(SETTINGS) $(APP_DATA)/config.yaml
	cp -r ./vendor/primer3_config $(APP_DATA) 
	cp -r ./assets/addgene/db/** $(APP_DATA) 
	cp -r ./assets/igem/db/** $(APP_DATA)
	cp -r ./assets/dnasu/db/** $(APP_DATA)
	cp ./assets/snapgene/features.tsv $(APP_DATA)
	cp ./assets/neb/enzymes.tsv $(APP_DATA)

ifeq ($(PLATFORM),Linux)
	install ./bin/linux $(APP)
	install -C ./vendor/linux/blastn $(LOCAL_BIN)
	install -C ./vendor/linux/ntthal $(LOCAL_BIN)
	install -C ./vendor/linux/primer3_core $(LOCAL_BIN)
	install -C ./vendor/linux/blastdbcmd $(LOCAL_BIN)
endif

ifeq ($(PLATFORM),Darwin)
	install ./bin/darwin $(APP)
	install -C ./vendor/darwin/blastn $(LOCAL_BIN)
	install -C ./vendor/darwin/ntthal $(LOCAL_BIN)
	install -C ./vendor/darwin/primer3_core $(LOCAL_BIN)
	install -C ./vendor/darwin/blastdbcmd $(LOCAL_BIN)
endif

windows: build
	cd scripts && makensis windows_installer.nsi

all: build install

dbs:
	cd assets && sh makeblastdbs.sh

uninstall: clean
	rm $(APP)
	rm -rf $(APP_DATA)

test:
	go test -timeout 200s ./internal/plade


dist: windows
	mkdir -p ${DIST}
	mkdir -p ${DIST_SRC}
	rsync -r --delete\
	 --exclude={'.git','dist','test','scripts','bin/plade_install.exe','bin/plade.exe','vendor/windows','assets/addgene/addgene.json','assets/dnasu/DNASU*','assets/igem/xml*','assets/neb/*/'}\
	 . ${DIST_SRC}
	tar -czf ${DIST_SRC_TAR} ${DIST_SRC}
	rm -rf ${DIST_SRC}

	zip ${DIST_WIN_ZIP} ./bin/plade_install.exe

	scp ${DIST_SRC_TAR} jjtimmons@frs.sourceforge.net:/home/frs/project/plasmid-defragger/
	scp ${DIST_WIN_ZIP} jjtimmons@frs.sourceforge.net:/home/frs/project/plasmid-defragger/