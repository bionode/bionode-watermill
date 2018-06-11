FROM bionode/bionode-watermill:dev

# tool installs taken from biocontainers

# === BOWTIE2 ===

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

# === BWA ===

RUN conda install bwa=0.7.15

# === SAMTOOLS ===

RUN conda install samtools=1.3.1

# === SRA-TOOLKIT ===

USER root
WORKDIR /bin/
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz
RUN tar -zxvf sratoolkit.2.8.2-1-ubuntu64.tar.gz
RUN rm -rf sratoolkit.2.8.2-1-ubuntu64.tar.gz
ENV PATH="${PATH}:/bin/sratoolkit.2.8.2-1-ubuntu64/bin"

# === bionode-ncbi ===

# TODO figure out how to do this another way
# TODO pre-accept the anonymous tracking
RUN chown root:biodocker /usr/local/lib/node_modules
RUN npm install -g bionode-ncbi
RUN chown -R biodocker:biodocker /home/biodocker/.config

# === FINISH ===

USER biodocker
