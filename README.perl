The perl based version of Methpup extends the classical version of Methpup as follows.

1. Full compute farm support with checkpointing, error logging, and file organisation.
2. Demultiplexing of reads.
3. Easier parameter input via ini files. Allows storage and replay of exact pipeline settings. 
4. Methplot graphs drawn as standard for each gene/well intersection detected.
5. Bisulphite conversion efficiency plots
6. Read drop out estimation at each step for each gene/well intersection detected.
7. Various output flat files to facillitate downstream analysis.

===================== 
Pre-requirement:
=====================

See dependencies.tab for particular version requirements of pipeline tools.
This version no longer requires python be installed or methsite.py and indel.py to be copied.

Additionally requires

R (tested with version 3.1.0) with libraries reshape and ggplot2 installed.


====================
Synopsis
====================

perl ./Methpup.pl -[i]ni_file -[h]elp


====================
Description
====================

The options for the program as as follows:

	-[i]ni_file		Path to ini file containing run configuration settings. A stub is available in ini/default.ini of this repository. 
	-[h]elp		This message.
	
====================
Test
===================

TODO


