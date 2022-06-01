FROM python:3.8

MAINTAINER "Pradeep Bashyal" pbashyal@nmdp.org

# Setup Repo for Neo4j version 3.5
RUN wget -O - https://debian.neo4j.com/neotechnology.gpg.key | apt-key add - && \
    echo 'deb https://debian.neo4j.com stable 3.5' | tee /etc/apt/sources.list.d/neo4j.list

# Install neo4j and cpanm
RUN apt-get update && \
  apt-get -y install neo4j cpanminus && \
  apt-get autoremove

# Install the Perl libraries
RUN cpanm install --auto-cleanup 0 --notest \
    REST::Client JSON JSON::Parse Math::Round

# Disable PIP Version check
ENV PIP_DISABLE_PIP_VERSION_CHECK=1
# Disable cache
ENV PIP_NO_CACHE_DIR=1

# Install Python libraries
RUN pip install pandas==1.0.3 networkx==2.3 py2neo==4.3.0

# Map cloned grimm directory to /grimm
# -v $PWD:/grimm
VOLUME /grimm
WORKDIR /grimm

# Setup env for Neo4J
ENV NEO4J_HOME=/usr/share/neo4j
ENV JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64

COPY docker-entrypoint.sh /usr/local/bin/

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

