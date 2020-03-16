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

COPY install_latex.R .

RUN Rscript install_latex.R

# Install other custom software

ENV APPS_HOME=/apps
RUN mkdir $APPS_HOME
WORKDIR $APPS_HOME

### gnparser ###
ENV APP_NAME=gnparser
ENV VERSION=0.13.1
ENV DEST=$APPS_HOME/$APP_NAME/$VERSION
RUN wget https://gitlab.com/gogna/gnparser/uploads/481ea1f2e32362fce661dc82b41cef36/$APP_NAME-v$VERSION-linux.tar.gz \
  && tar xf $APP_NAME-v$VERSION-linux.tar.gz \
  && rm $APP_NAME-v$VERSION-linux.tar.gz \
  && mv "$APP_NAME" /usr/local/bin/
  
WORKDIR /home/rstudio
