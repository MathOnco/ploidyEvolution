# Conda Environment For HATCHet

### Contents 

1. [Requirements](#requirements)
1. [HATCHet Documentation](#hatchet-documentation)
1. [Gurobi Version](#gurobi-version)
1. [Gurobi License](#gurobi-license)
1. [Build Docker Image](#build-docker-image)
1. [Start Interactive Docker Container](#start-interactive-docker-container)
1. [Create Singularity Recipe from Dockerfile](#create-singularity-recipe-from-dockerfile)
1. [Build Singularity Image for HPC](#build-singularity-image-for-hpc)
1. [Start Interactive Singularity Container](#start-interactive-singularity-container)
1. [Execute Scripts in Container](#execute-scripts-in-container)
1. [Running Scripts in Docker Container Non Interactive](#running-scripts-in-docker-container-non-interactive)
1. [Running Scripts in Singularity Container Non Interactive](#running-scripts-in-singularity-container-non-interactive)
1. [Example Shell Script Input](#example-shell-script-input)

### Requirements

1. Docker (built with v2.2.0.5)
1. Gurobi v9.0.2
1. Python 3+

### HATCHet Documentation

1. [Docs](https://github.com/raphael-group/hatchet)

### Gurobi Version

1. This project was created with Gurobi v9.0.2
1. **IMPORTANT**: If another version of gurobi is used be sure to update the Dockerfile and Singularity recipe to reflect the correct version number
1. **IMPORTANT** The Gurobi .tar.gz file **must** be placed in **hatchetContainers/gurobi** in order to build the Docker/Singularity images successfully

### Gurobi License

1. Gurobi is part of the HATCHet pipeline and requires a license to use.
1. To include your license in the Docker/Singularity containers you must mount the directory containing the license, as described in running the containers in the following sections
1. The gurobi license must be named **gurobi.lic** (as recommended by Gurobi)
1. If mounted correctly the license can be accessed inside of the containers by the environment variable **GRB_LICENSE_FILE** 

Example:
```sh 
$ cat $GRB_LICENSE_FILE
# OUTPUT - you should see the contents of gurobi.lic
```

### Build Docker Image

1. **IMPORTANT**: Do not forget the period at the end of the docker build command

```sh 
$ cd hatchetContainers
$ docker build -t hatchet:<VERSION_NUMBER> .
```
Example:
```sh 
$ docker build -t hatchet:1.0 .
```

### Start Interactive Docker Container

```sh 
$ docker run -it --rm -v /path/to/local/input:/home/docker/input -v /path/to/local/output:/home/docker/output -v /path/to/local/license:/home/docker/license hatchet:<VERSION_NUMBER>
```
1. The Docker container contains two directories that can be mapped to local directories for mount points
	- **-it**: allows the container to be attached an run interactively
	- **--rm**: shutdown and remove container after execution
	- **-v**: allows the user to set these mount points 
1. **IMPORTANT**: The internal mount points are */home/docker/input*, */home/docker/output* and */home/docker/license*.  These cannot be changed without editing the Dockerfile and rebuilding the image

### Create Singularity Recipe from Dockerfile

```sh 
$ cd hatchetContainers
$ python3 docker2singularity.py 
```

### Build Singularity Image for HPC

1. **IMPORTANT**: To build a Singularity image root access is required as well as a Linux distro.  This cannot be done on the HPC, but can be on a local instance of Linux.
```sh
$ cd hatchetContainers
$ sudo singularity build <IMAGE_NAME>.sif <SINGULARITY_RECIPE>
```
1. Example:
```sh
$ sudo singularity build hatchet.sif Singularity
```
1. You may now copy this **.sif** image file to the HPC

### Start Interactive Singularity Container

1. **NOTE**: The Singularity container will have access to your home directory automatically
1. **NOTE**: Since space is limited in your home directory it is advised to mount directories to the input and output directories provided in the container using the **-B** flag, as shown above
```sh 
$ module load singularity/3.10
$ singularity shell -B /path/to/local/input:/home/docker/input -B /path/to/local/output:/home/docker/output -B /path/to/local/license:/home/docker/license hatchet.sif 
```
1. Once the container is activated run these commands to start the built-in Conda environment:
```sh 
$ . "/opt/conda/etc/profile.d/conda.sh"
$ conda activate hatchet
```

### Execute Scripts in Container

1. Executing scripts inside of either a Docker container or a Singularity Container is identical
1. Scripts can be run as if running them on your local machine

### Running Scripts in Docker Container Non Interactive

```sh 
$ docker run --rm -v /path/to/local/input:/home/docker/input -v /path/to/local/output:/home/docker/output -v /path/to/local/license:/home/docker/license hatchet:<VERSION_NUMBER> bash <SHELL_SCRIPT>
```
1. Parameters:
	- **--rm**: shutdown and remove container after execution
	- **-v**: allows the user to set these mount points 
1. Example:
```sh 
docker run --rm -v /path/to/local/input:/home/docker/input -v /path/to/local/output:/home/docker/output -v /path/to/local/license:/home/docker/license hatchet:1.0 bash /home/docker/input/example.sh
```
1. **NOTE**: See below for [example.sh](#example-shell-script-input) content

### Running Scripts in Singularity Container Non Interactive

1. If running Sinuglarity on the HPC first load the Singularity module:
```sh 
$ module load singularity/3.10
```

```sh 
$ singularity exec -B path/to/local/input:/home/docker/input -B /path/to/local/output:/home/docker/output -B /path/to/local/license:/home/docker/license <SINGULARITY_CONTAINER>.sif bash <SHELL_SCRIPT>
```
1. Parameters:
	- **-B**: allows the user to set these mount points 
1. Example:
```sh 
$ singularity exec -B /share/NGS -B /path/to/local/input:/home/docker/input -B /path/to/local/output:/home/docker/output -B /path/to/local/license:/home/docker/license hatchet.sif bash example.sh
```
1. **NOTE**: See below for [example.sh](#example-shell-script-input) content

### Example Shell Script Input

1.  It is important to add these two lines to the top of each of your shell scripts so the container can source the internal Conda Environment
```text 
. "/opt/conda/etc/profile.d/conda.sh"
conda activate hatchet
``` 

```text 
#!/bin/bash
# example.sh

. "/opt/conda/etc/profile.d/conda.sh"
conda activate hatchet

min=0
max=5

python /opt/hatchet/utils/BBot.py \
	/home/docker/input/P06269.tune.bbc \
	--segfile /home/docker/input/P06269.tune.seg \
	--command RD \
	--ymin $min \
	--ymax $max \
	-x /home/docker/output
```
