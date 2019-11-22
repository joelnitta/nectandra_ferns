# Only run this after making renv.lock by
# running install_packages.R

FROM rocker/verse:3.6.1

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
  mafft \
  iqtree \
  fasttree

# Install R packages using Renv snapshot

COPY ./renv.lock ./

COPY ./renv_restore.R ./

RUN mkdir renv

RUN Rscript renv_restore.R

RUN echo '.libPaths("/renv")' >> /usr/local/lib/R/etc/Rprofile.site

# Install latex packages with tinytex

RUN Rscript install_latex.R

# Install other custom software

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### gnparser ###
ENV APP_NAME=gnparser
ENV VERSION=0.7.5
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://www.dropbox.com/s/7jcrjj0o39vuh3x/$APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && rm $APP_NAME-v$VERSION-linux.tar.gz?dl=1 \
  && mv "$APP_NAME" /usr/local/bin/
