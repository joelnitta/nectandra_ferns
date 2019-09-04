# Building the image

This is how to build the docker image.

First, launch the `rocker:verse` container, and run `install_packages.R`. This will install packages to a private library with `renv` and (more importantly) write the `renv.lock` file, which contains version information for all packages installed.
  
```
docker run --rm -e DISABLE_AUTH=true -v /Users/joelnitta/repos/nectandra:/home/rstudio/project rocker/verse:3.6.1 bash
Rscript install_packages.R
exit
```

Once that's done, build and tag the image.

```
docker build . -t joelnitta/nectandra_ferns:3.6.1
```
