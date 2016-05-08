#!/bin/sh
rm MatlabSignalServer* 
rm *.html
cd /afs/fusione.it/user/e/eneatesi/Daniele_Carnevale/IRE_REF

MARTESERVER=boncagni:8040
sleep 10
wget --user-agent=Mozilla/5.0  http://$MARTESERVER/BROWSE/StateMachine/?StatusChangeRequest=PULSE_SETUP_COMPLETED
sleep 1
wget --user-agent=Mozilla/5.0  http://$MARTESERVER/BROWSE/StateMachine/?StatusChangeRequest=PRE
sleep 20
wget --user-agent=Mozilla/5.0  http://$MARTESERVER/BROWSE/StateMachine/?StatusChangeRequest=EJP
sleep 1
wget --user-agent=Mozilla/5.0  http://$MARTESERVER/BROWSE/StateMachine/?StatusChangeRequest=COLLECTION_COMPLETED
sleep 3 
wget --user-agent=Mozilla/5.0  --post-data 'ALLSIGNALS=Yes&FORMSENT=Yes' http://$MARTESERVER/BROWSE/MatlabSignalServer
mv MatlabSignalServer signalsremote.mat
rm index.html*
sync
