# Docker 
In order to use the dockerfiles and/or docker containers, Docker Desktop must be installed. Docker [installation](https://www.docker.com/products/docker-desktop/) is straightforward. 
Once installed, Docker engine can be started using the graphic interface of Docker. When the engine is running, the containers can be built and run using the terminal's commands.
It is suggested to change the resources that Docker engine is using, since containers may require more memory than the default allocation. In order to do this, open Docker Desktop, then click on settings->resources, adjust based on personal needs. It is also useful to change the default storage limite of the Docker Engine in settings->Docker Engine and change the builder storage default.

## Docker containers
A container is an isolated environment based on a specific operating system (OS) (typically Linux based), in which a solver can be built. The environment created in a container includes all the dependencies of the solver that needs to be built. Once the container is created, it can be run on any machine and the solver can be built directly in this virtual environment. 
In this work the containers are created using software [Docker](https://www.docker.com). A Docker container is built from a dockerfile. After it is created, it is saved as an image that can be copied anywhere is needed or uploaded on the cloud [DockerHub](https://hub.docker.com).

## Dockerfile
In the folder Docker/ there are three subfolders solver/, ubuntu20/, ubuntu22/ containing a dockerfile each. A dockerfile usually begins with an image import of the OS, which in this case is Ubuntu-20.04 and Ubuntu-22.04. It could also start by importing an image that contain the desired environment such as in the dockerfile in the solver/ folder.
For more details about the dockerfiles created in this work, refer to the [simvascular DockerHub page](https://registry.hub.docker.com/u/simvascular).

## Build a container
In this section the steps to build the image with a pre-compiled svFSIplus solver are briefly described. To create the image from the dockerfile provided in Docker/solver, follow the steps below:
1) build an Ubuntu-based image containing the whole environment in which svFSIplus program can be compiled. The provided dockerfiles are based on Ubuntu-20.04 and Ubuntu-22.04, but they can be easily adapted to use the latest version of Ubuntu, by changing the following line in Docker/ubuntu20/dockerfile or Docker/ubuntu22/dockerfile: 
```
FROM ubuntu:20.04 AS base / FROM ubuntu:22.04 AS base
```
to 
```
FROM ubuntu:latest AS base
```
Build the environmnet Docker image: 
```
cd Docker/ubuntu20 or cd Docker/ubuntu22
```
```
docker build -t RepositoryName:tagImage .
```
where -t allows the user to set the name for the image created. For example:
```
docker build -t libraries:latest .
```
2) build the image containing the compiled svFSIplus program. This image will be based on the environment created in the previous step (libraries:latest). In order to do this, open the Docker/solver/dockerfile and modify the following lines:
```
FROM simvascular/libraries:ubuntu22 AS builder 
```
to
```
FROM libraries:latest AS builder
```
and 
```
FROM simvascular/libraries:ubuntu22 AS final 
```
to
```
FROM libraries:latest AS final
```
Build the solver image:
```
cd Docker/solver
```
```
docker build -t solver:latest .
```
The image include the PETSc-based svFSIplus executable in:
```
/build-petsc/svFSIplus-build/bin/svfsiplus
```
and the Trilinos-based svFSIplus executable in:
```
/build-trilinos/svFSIplus-build/bin/svfsiplus
```
## Run a container
Once the image is created, it can be run interactively by running the following command:

docker run -it -v FolderToUpload:/NameOfFolder RepositoryName:tagImage

In this command:
- -it: means run interactively Docker image
- -v: mounts a directory 'FolderToUpload' from the host machine in the container where the directory has the name '/NameOfFolder'. For example the directory containing the source code of the solver and all the test cases can be mounted in the container.

When the container is running interactively, the solver can be built and the test cases can be run. When exiting the container, all the output files inside '/NameOfFolder' are not deleted and they will still be available in the host machine.
