FROM python:3.9.13

MAINTAINER Cecilia Jensen "cecilia.jensen@tum.de"

LABEL website=https://gitlab.lrz.de/proteomics/topas/wp3_sample_pipeline
LABEL description="WP3 sample pipeline"
LABEL tags="mass spectrometry tmt isobaric labeling"
LABEL documentation=https://gitlab.lrz.de/proteomics/topas/wp3_sample_pipeline

# Tell docker that we don't want to be bothered with questions
ARG DEBIAN_FRONTEND=noninteractive

# set root directory
ENV HOME /root
WORKDIR /root

# for mono installation
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF
RUN echo "deb https://download.mono-project.com/repo/debian stable-buster main" | tee /etc/apt/sources.list.d/mono-official-stable.list

RUN apt-get update && apt-get install -y \
        mono-devel \
        zip \
    && rm -rf /var/lib/apt/lists/*

# set root directory
ENV HOME /root
WORKDIR /root

RUN pip install poetry==1.5.1
# poetry uses virtualenvs by default -> we want global installation
RUN poetry config virtualenvs.create false
ADD pyproject.toml /root/pyproject.toml
ADD poetry.lock /root/poetry.lock
RUN poetry install --no-cache

COPY . /root/wp3_sample_pipeline
ADD hash.file /root/hash.file

WORKDIR /root/wp3_sample_pipeline

