#!/bin/bash

# Make sure to add the following entry in /etc/default/grub:
#GRUB_CMDLINE_LINUX="cgroup_enable=memory swapaccount=1"



# Name
NAME="korseby/metfamily-dev"

# CPU options
#CPU_SHARES="--cpu-shares=8"
#CPU_SETS="--cpuset-cpus=0-$[$CPU_SHARES-1]"
#CPU_MEMS="--cpuset-mems=0"
#MEM="--memory=8g"

# Ports
PORT_PUB=9011
PORT_DOCKER=3838

# Volumes
VOL="--volume=/vol/R/shiny/srv/shiny-server/MetFam:/vol/R/shiny/srv/shiny-server/MetFam:ro"



# Run docker
docker run --publish=${PORT_PUB}:${PORT_DOCKER} --log-driver=syslog $VOL $CPU_SHARES $CPU_SETS $CPU_MEMS $MEM --name="$(echo ${NAME} | sed -e 's/.*\///')-run" -i -t -d $NAME



# Detach/Attach docker
# detach: CTRL-P + CTRL-Q

# Start and attach docker (you can also use docker start -ai instead)
#docker start ${NAME}-run
#docker attach ${NAME}-run

# Start shell inside running docker
#docker exec -i -t ${NAME}-run /bin/bash

# Start failed container with different entrypoint
#docker run -ti --entrypoint=/bin/bash ${NAME}-run

# Commit changes locally
#docker commit ${NAME}-run

# Show docker container and images
#docker ps -a
#docker images

# Delete container and image
#docker rm ${NAME}-run
#docker rmi ${NAME}

# Delete exited containers
#docker rm $(docker ps -a -f status=exited | grep -v CONTAINER\ ID | awk '{print $1}')

# Delete intermediate/untagged images
#docker rmi $(docker images | grep '^<none>' | awk '{print $3}')

