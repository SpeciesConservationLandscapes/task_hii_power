FROM scl3/task_base:latest

RUN pip install git+https://github.com/SpeciesConservationLandscapes/task_base.git \
    && pip install pandas==1.3.2 \
    && pip install scikit-learn==1.0.1

WORKDIR /app
COPY $PWD/src .
