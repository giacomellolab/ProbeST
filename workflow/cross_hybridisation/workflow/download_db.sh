#!/bin/bash
# List of files to download
files=(

env_nt-nucl-metadata.json
env_nt.tar.gz
env_nt.tar.gz.md5
nt_prok-nucl-metadata.json
nt_prok.00.tar.gz
nt_prok.00.tar.gz.md5
nt_prok.01.tar.gz
nt_prok.01.tar.gz.md5
nt_prok.02.tar.gz
nt_prok.02.tar.gz.md5
nt_prok.03.tar.gz
nt_prok.03.tar.gz.md5
nt_prok.04.tar.gz
nt_prok.04.tar.gz.md5
nt_prok.05.tar.gz
nt_prok.05.tar.gz.md5
nt_prok.06.tar.gz
nt_prok.06.tar.gz.md5
nt_prok.07.tar.gz
nt_prok.07.tar.gz.md5
nt_prok.08.tar.gz
nt_prok.08.tar.gz.md5
nt_prok.09.tar.gz
nt_prok.09.tar.gz.md5
nt_prok.10.tar.gz
nt_prok.10.tar.gz.md5
nt_prok.11.tar.gz
nt_prok.11.tar.gz.md5
nt_prok.12.tar.gz
nt_prok.12.tar.gz.md5
nt_prok.13.tar.gz
nt_prok.13.tar.gz.md5
nt_prok.14.tar.gz
nt_prok.14.tar.gz.md5
nt_prok.15.tar.gz
nt_prok.15.tar.gz.md5
nt_prok.16.tar.gz
nt_prok.16.tar.gz.md5
nt_prok.17.tar.gz
nt_prok.17.tar.gz.md5
nt_prok.18.tar.gz
nt_prok.18.tar.gz.md5
nt_prok.19.tar.gz
nt_prok.19.tar.gz.md5
nt_prok.20.tar.gz
nt_prok.20.tar.gz.md5
nt_viruses-nucl-metadata.json
nt_viruses.00.tar.gz
nt_viruses.00.tar.gz.md5
nt_viruses.01.tar.gz
nt_viruses.01.tar.gz.md5
nt_viruses.02.tar.gz
nt_viruses.02.tar.gz.md5
nt_viruses.03.tar.gz
nt_viruses.03.tar.gz.md5
nt_viruses.04.tar.gz
nt_viruses.04.tar.gz.md5
nt_viruses.05.tar.gz
nt_viruses.05.tar.gz.md5
nt_viruses.06.tar.gz
nt_viruses.06.tar.gz.md5
nt_viruses.07.tar.gz
nt_viruses.07.tar.gz.md5
nt_viruses.08.tar.gz
nt_viruses.08.tar.gz.md5
nt_viruses.09.tar.gz
nt_viruses.09.tar.gz.md5
nt_viruses.10.tar.gz
nt_viruses.10.tar.gz.md5
nt_viruses.11.tar.gz
nt_viruses.11.tar.gz.md5
nt_viruses.12.tar.gz
nt_viruses.12.tar.gz.md5
nt_viruses.13.tar.gz
nt_viruses.13.tar.gz.md5
nt_viruses.14.tar.gz
nt_viruses.14.tar.gz.md5
nt_viruses.15.tar.gz
nt_viruses.15.tar.gz.md5
nt_viruses.16.tar.gz
nt_viruses.16.tar.gz.md5
nt_viruses.17.tar.gz
nt_viruses.17.tar.gz.md5
nt_viruses.18.tar.gz
nt_viruses.18.tar.gz.md5
nt_viruses.19.tar.gz
nt_viruses.19.tar.gz.md5
nt_viruses.20.tar.gz
nt_viruses.20.tar.gz.md5 

)

# FTP base URL
base_url="https://ftp.ncbi.nlm.nih.gov/blast/db/"

# Loop through each file and download it
for file in "${files[@]}"; do
    wget "${base_url}${file}"
done

