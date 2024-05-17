# Docker 
In order to use the dockerfiles and/or docker containers, Docker Desktop must be installed. Docker [installation](https://www.docker.com/products/docker-desktop/) is straightforward. 
Once installed, Docker engine can be started using the graphic interface of Docker. When the engine is running, the containers can be built and run using the terminal's commands.
It is suggested to change the resources that Docker engine is using, since containers may require more memory than the default allocation. In order to do this, open Docker Desktop, then click on settings->resources, adjust based on personal needs. It is also useful to change the default storage limite of the Docker Engine in settings->Docker Engine and change the builder storage default.

## Docker containers
A container is an isolated environment based on a specific operating system (OS) (typically Linux based), in which a solver can be built. The environment created in a container includes all the dependencies of the solver that needs to be built. Once the container is created, it can be run on any machine and the solver can be built directly in this virtual environment. 
<<<<<<< HEAD
In this work the containers are created using [Docker](https://www.docker.com). A Docker container is built from a dockerfile. After it is created, it is saved as an image that can be copied anywhere is needed or uploaded on the cloud [DockerHub](https://hub.docker.com).
=======
In this work the containers are created using software [Docker](https://www.docker.com). A Docker container is built from a dockerfile. After it is created, it is saved as an image that can be copied anywhere is needed or uploaded on the cloud [DockerHub](https://hub.docker.com).
>>>>>>> b887fd6aa14382cd941909a1ce9fedcf985c2ed3

## Dockerfile
In the folder Docker/ there are two subfolders ubuntu20/ and ubuntu22/ containing a dockerfile each. A dockerfile usually begins with an image import of the environment OS, which in this case is Ubuntu-20.04 and Ubuntu-22.04. 
For more details about the dockerfiles created in this work, refer to the DockerHub page [personalDockerHub](https://hub.docker.com/repository/docker/dcodoni/lib/general).

## Build a container
Once the dockerfile is created, it must be built to oget the final image of the container. This can be done by moving to the directory where dockerfile is and run the command:

docker buildx build -t RepositoryName:tagImage .

where -t allows the user to set the name for the image created. 
If your file is not named dockerfile, then in order to build a particular file just use:

docker buildx build -f filename -t RepositoryName:tagImage .

## Run a container
Once the image is created, it can be run interactively by running the following command:

docker run -it -v FolderToUpload:/NameOfFolder RepositoryName:tagImage

In this command:
- -it: means run interactively Docker image
- -v: mounts a directory 'FolderToUpload' from the host machine in the container where the directory has the name '/NameOfFolder'. For example the directory containing the source code of the solver and all the test cases can be mounted in the container.

When the container is running interactively, the solver can be built and the test cases can be run. When exiting the container, all the output files inside '/NameOfFolder' are not deleted and they will still be available in the host machine.
