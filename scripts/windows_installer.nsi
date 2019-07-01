!include EnvVarUpdate.nsh

; installer's name
Name "repp Installer"

; name of the install binary
OutFile "..\bin\repp_install.exe"

; installation directory
InstallDir $PROFILE\.repp

; create an installer exe
Section "install"

; add the directory to Path
${EnvVarUpdate} $0 "PATH"  "A" "HKCU" $INSTDIR

; set the output path to the installation directory
SetOutPath $INSTDIR

; the binary
File "..\bin\repp.exe"

; configuration
File "..\config\config.yaml"

; primer3 configuration directory
SetOutPath "$INSTDIR\primer3_config"
File /nonfatal /r "..\vendor\primer3_config\"
SetOutPath $INSTDIR

; BLAST dbs
File /nonfatal /r "..\assets\addgene\db\"
File /nonfatal /r "..\assets\igem\db\"
File /nonfatal /r "..\assets\dnasu\db\"

; feature and enzyme dbs
File "..\assets\snapgene\features.tsv"
File "..\assets\neb\enzymes.tsv"

; blastn, primer3_core, ntthal, other dlls
File /nonfatal /r "..\vendor\windows\"

SectionEnd

; create an uninstaller exe
Section "uninstall"

; remove the directory from the user's Path
${un.EnvVarUpdate} $0 "PATH"  "R" "HKCU" $INSTDIR

SectionEnd