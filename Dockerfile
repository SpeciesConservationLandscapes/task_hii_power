FROM scl3/task_base:latest

RUN pip install git+https://github.com/SpeciesConservationLandscapes/task_base.git

WORKDIR /app
COPY $PWD/src .
