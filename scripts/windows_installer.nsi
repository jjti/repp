!include EnvVarUpdate.nsh

; installer's name
Name "rvec Installer"

; name of the install binary
OutFile "..\bin\rvec_install.exe"

; installation directory
InstallDir $PROFILE\.rvec

; stuff to install
Section ""

; add the directory to Path
${EnvVarUpdate} $0 "PATH"  "A" "HKCU" $INSTDIR

; set the output path to the installation directory
SetOutPath $INSTDIR

; the binary
File "..\bin\rvec.exe"

; configuration
File "..\config\config.yaml"

; primer3 configuration directory
SetOutPath "$INSTDIR\primer3_config"
File /nonfatal /a /r "..\vendor\primer3_config\"
SetOutPath $INSTDIR

; BLAST dbs
File /nonfatal /a /r "..\assets\addgene\db\"
File /nonfatal /a /r "..\assets\igem\db\"
File /nonfatal /a /r "..\assets\dnasu\db\"

; feature and enzyme dbs
File "..\assets\snapgene\features.tsv"
File "..\assets\neb\enzymes.tsv"

; blastn, primer3_core, ntthal, other dlls
File /nonfatal /a /r "..\vendor\windows\"

SectionEnd