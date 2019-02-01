#!/bin/bash

# Name
NAME="korseby/metfamily-dev"

# Cache
#CACHE="--no-cache"

# CPU options
#CPU_SHARES="--cpu-shares=8"
#CPU_SETS="--cpuset-cpus=0-$[$CPU_SHARES-1]"
#CPU_MEMS="--cpuset-mems=0"
#MEM="--memory=8g"



# Build docker
docker build $CACHE --rm=true $CPU_SHARES $CPU_SETS $CPU_MEMS $MEM --tag=$NAME .

#docker push $NAME

