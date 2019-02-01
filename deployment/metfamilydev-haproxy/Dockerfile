FROM debian:jessie

# Install haproxy
RUN sed -i 's/^# \(.*-backports\s\)/\1/g' /etc/apt/sources.list
RUN apt-get -y update
RUN apt-get -y install haproxy
RUN sed -i 's/^ENABLED=.*/ENABLED=1/' /etc/default/haproxy

# Add configuration files
ADD haproxy.cfg /etc/haproxy/haproxy.cfg
ADD start.sh /haproxy-start

# Define volumes
VOLUME ["/haproxy-override"]

# Expose ports
EXPOSE 80
EXPOSE 443

# Working directory
WORKDIR /etc/haproxy

# Entrypoint
ENTRYPOINT ["/bin/bash", "/haproxy-start"]

