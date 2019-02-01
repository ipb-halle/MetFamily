#!/bin/sh

# start
docker-compose up -d
docker-compose scale metfamilydev=8

# status
#docker-compose ps
#docker-compose logs

# stop
#docker-compose down

