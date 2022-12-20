FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
MAINTAINER Chris Plaisier <plaisier@asu.edu>
RUN apt-get update

RUN apt-get install --yes \
 vim-common \
 python3 \
 python3-pip \
 git

RUN pip3 install pandas numpy tdqm statsmodels

RUN git clone https://github.com/plaisier-lab/OncoMerge.git

