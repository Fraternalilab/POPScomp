# History of Uploads
The plan is to perform automated weekly updates.
The PDB is updated weekly, Wednesday at 00:00 UTC,
therefore the POPS update can be run on Wednesday morning,
for example at 3:00.

## 
Increment

## Wed 15 Jan 23:07:43 GMT 2020
Complete upload of JSON output for PDBML (<4MB).

```
jkleinj@ac:~/database$ lftp -e "mirror -R /home/jkleinj/database/JSONVAL/ /upload/" -u pops,2n8BDJiY ftp-private.ebi.ac.uk
Total: 1060 directories, 155575 files, 0 symlinks
New: 155575 files, 0 symlinks
52372354574 bytes transferred in 45738 seconds (1.09 MiB/s)
lftp pops@ftp-private.ebi.ac.uk:/>
```

