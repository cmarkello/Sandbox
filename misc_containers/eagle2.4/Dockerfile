FROM ubuntu:20.10

WORKDIR /usr/src/app

RUN apt-get update && apt-get install -y libgomp1 tabix

COPY Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz /usr/src/app/genetic_map_hg19_withX.txt.gz
COPY Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz /usr/src/app/genetic_map_hg38_withX.txt.gz
COPY Eagle_v2.4.1/eagle /usr/src/app/eagle
COPY Eagle_v2.4.1/src /usr/src/app/src

CMD [ "/bin/bash" ]

