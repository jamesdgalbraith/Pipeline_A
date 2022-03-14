# run RepeatModeler
### create database
BuildDatabase -name seq/${GENOME} seq/${GENOME}

### run repeatmodeler
RepeatModeler -pa ${THREADS} -database seq/${GENOME}

### remove database
rm seq/${GENOME}.n* seq/${GENOME}.translation