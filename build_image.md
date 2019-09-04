# Building the image

This is how to build the docker image.

First, launch the `rocker:verse` container: 
  
```
docker run --rm -e -d DISABLE_AUTH=true -v /Users/joelnitta/repos/nectandra:/home/rstudio/project rocker/verse:3.6.1
```

Get the name of container using `docker ps`, then enter the container and run `install_packages.R`. You may have to run `renv::consent()` first.

This will install packages to a private library with `renv` and (more importantly) write the `renv.lock` file, which contains version information for all packages installed.

```
docker exec -it <NAME OF CONTAINER> bash
cd /home/rstudio/project
Rscript install_packages.R
```

Once that's done, build and tag the image.

```
docker build . -t joelnitta/nectandra_ferns:3.6.1
```
