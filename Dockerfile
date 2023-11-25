# source image
FROM python:3.9

# set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
  apt-get update && \
  apt-get upgrade -y && \
  apt-get install -y  \
      build-essential \
      curl

# apt clean and remove cached source lists
RUN apt-get clean && \
  rm -rf /var/lib/apt/lists/*

# install poetry
RUN curl -sSL https://install.python-poetry.org | python3 -
ENV PATH="/root/.local/bin:$PATH"

# copy app code
COPY . /app
WORKDIR /app

# install python requirements
RUN poetry config virtualenvs.create false
RUN poetry install

# default command
CMD ["phenopy"]
