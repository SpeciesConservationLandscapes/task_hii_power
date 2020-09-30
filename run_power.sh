#!/bin/bash

docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Afrotropic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Australasia'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Indomalayan'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Nearctic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Neotropic'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Oceania'
docker run -it -v $PWD/.config:/root/.config scl3/task_hii_power python hii_power.py -r 'Palearctic'
