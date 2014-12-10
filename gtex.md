
### Get expression/genotype data from GTEx exchange area
1. Download from GTEx Exchange Area (Provisinal tab)
- Install aspera
- download from http://downloads.asperasoft.com/downloads (select aspera connect and then show all installers)
- Download using aspera
`> "/home/him1/.aspera/connect/bin/ascp" -QTr -l 300M -k 1 -i "/home/him1/.aspera/connect/etc/asperaweb_id_dsa.openssh" -W THIS_IS_COPIED_FROM_DBGAP dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/haeky2/41316 .`
`> "/userhome/anuar/.aspera/connect/bin/ascp" -QTr -l 300M -k 1 -i
/userhome/anuar/.aspera/connect/etc/asperaweb_id_dsa.openssh
/userhome/anuar/.aspera/connect/etc/asperaweb_id_dsa.openssh -W A95BE876F66774D76A47EB0C286B1EC2C0834A7F36BF5466E89703EE38E5088DB66BE17A6FB06466792BD90872536274ED dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/haeky2/41400 .`
2. Download manifest
3. Get repository key
4. Decrypt
- download decryption tools from dbGaP
- run sratoolkit.jar (> java -jar sratoolkit.jar)
- import repository key within sratookit.jar
- from the working directory run vbd-decrypt (> /group/im-lab/nas40t2/haky/Data/decryption.2.4.2-ubuntu64/bin/vdb-decrypt 41400/)
5. Note: had to download on genegate, rsync to tarbell and decrypt on tarbell
`> /group/im-lab/nas40t2/haky/Data/decryption.2.4.2-ubuntu64/bin/vdb-decrypt   41400/
`
