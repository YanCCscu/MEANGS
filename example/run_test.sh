wget -O SRR039541.3_1.clean.fq.gz https://ndownloader.figshare.com/files/22431854
wget -O SRR039541.3_2.clean.fq.gz https://ndownloader.figshare.com/files/22431857
../meangs.py -1 SRR039541.3_1.clean.fq.gz -2 SRR039541.3_2.clean.fq.gz -o HumanMito -t 16 -n 2000000 -i 300 --deepin
