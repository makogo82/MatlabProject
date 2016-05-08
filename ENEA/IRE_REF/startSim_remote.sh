#!/bin/tcsh

eval `ssh-agent -c`

ssh-add .ssh/id_rsa
cd /afs/fusione.it/user/e/eneatesi/Daniele_Carnevale/IRE_REF/matlabFileSource
scp  -i .ssh/id_rsa.pub  IreRef.h   luca@boncagni:/home/luca/code/MARTe/GAMs/ControllerGAM/IreRef.h
scp  -i .ssh/id_rsa.pub  IreRef.cpp luca@boncagni:/home/luca/code/MARTe/GAMs/ControllerGAM/IreRef.cpp
ssh  -i .ssh/id_rsa.pub luca@boncagni  " cd /afs/fusione.it/user/e/eneatesi/code/MARTe/GAMs/ControllerGAM; rm ControllerGAMClassInfo.sinfo.cpp"
scp  -i .ssh/id_rsa.pub  /afs/fusione.it/user/e/eneatesi/code/MARTe/GAMs/ControllerGAM/ControllerGAMClassInfo.sinfo.cpp luca@boncagni:/home/luca/code/MARTe/GAMs/ControllerGAM/ControllerGAMClassInfo.sinfo.cpp

cd /afs/fusione.it/user/e/eneatesi/Daniele_Carnevale/IRE_REF
scp  -i .ssh/id_rsa.pub   data_actual.dat luca@boncagni:/home/luca/code/MARTe/TestArea/data_actual.dat
scp  -i .ssh/id_rsa.pub   zero_actual.dat luca@boncagni:/home/luca/code/MARTe/TestArea/zero_actual.dat
scp  -i .ssh/id_rsa.pub $1 luca@boncagni:/home/luca/code
ssh  -i .ssh/id_rsa.pub luca@boncagni  "cd /home/luca/code/MARTe/GAMs/ControllerGAM/; make -f Makefile.linux"
sleep 1
ssh  -i .ssh/id_rsa.pub luca@boncagni  "sync;killall MARTe.ex; cd /home/luca/code/MARTe/TestArea;  ./MARTe-frascati.sh ../../$1" &
./download_remote.sh 
ssh  -i .ssh/id_rsa.pub luca@boncagni  "killall MARTe.ex;"
 