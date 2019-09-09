# source image
FROM python:3.7

# set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y  \
      build-essential

# apt clean and remove cached source lists
RUN apt-get clean && \
  rm -rf /var/lib/apt/lists/*

# install pipenv
RUN pip install pipenv --upgrade

# copy app code
COPY . /app
WORKDIR /app

# install python requirements
RUN pipenv install --dev --deploy --system

# install phenopy
RUN pip install .

# default command
CMD ["phenopy"]
