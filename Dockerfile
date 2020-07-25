FROM ubuntu
RUN apt-get update
RUN apt-get install -y libtbb-dev
RUN mkdir /opt/accel-align-release/
WORKDIR /opt/accel-align-release/
COPY accalign-x86-64 .
COPY accindex-x86-64 .
