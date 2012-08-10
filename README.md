MatlabLiveWire
==============

LiveWire on matlab. 
Using IvusSnakes (http://ivussnakes.sourceforge.net/) ImageJ plugin to do the calculations. 
The initial version was a port of the aforementioned plugin, but it was too slow for any useful application (removed at f33a5a6 commit) and I went back to using Java to do the hard lifting.

Ended rewriting the edge costs and thus effectively implementing William A. Barrett* and Eric N. Mortensen Interactive live-wire boundary extraction. Medical Image Analysis (1996/7) volume 1, number 4, pp 331–341.

Written by Timo Rantalainen tjrantal(at)gmail.com 2012


Brief how to:

Obtaining the Code for Windows users (the last two steps are the same for linux users as well..):
Install git (http://git-scm.com/downloads/), just the least invasive installation suffices.
Open Git Bash
cd to the folder you'd like to copy the code to (Windows N.B. if you need to change hard drive you cd /c/folder, where c is your hard drive.)
git init													(which initializes the folder you are currently in)
git clone https://github.com/tjrantal/MatlabLiveWire.git	(copies the code from github)

Compiling the code:
open command line (cmd.exe) and cd to the folder where you git cloned in
javac javaEngineLiveWire/*.java

For testing the plugin, take a look at the testJava.m to load a DICOM image to test the liveWire on.
