import os
from spython.main.parse.parsers import DockerParser
from spython.main.parse.writers import SingularityWriter

TEMP_DOCKERFILE = 'tDockerfile'

###############################################################################
###############################################################################
# functions

def prep_dockerfile ():
    # convert ARG to ENV; ARG is not supported by Singularity
    tDockerfile = []
    passthru = [ tDockerfile.append(row.strip('\n').replace('ARG ', 'ENV ')) for row in open('./Dockerfile') ]

    with open(TEMP_DOCKERFILE, 'w') as f:
        for row in tDockerfile:
            f.write(row)
            f.write('\n')

def parse_and_write_singularity ():
    parser = DockerParser(TEMP_DOCKERFILE)
    writer = SingularityWriter(parser.recipe)
    
    with open('./Singularity', 'w') as f:
        f.write(writer.convert())

###############################################################################
###############################################################################
###############################################################################
# main()

if __name__ == '__main__':

    prep_dockerfile()
    parse_and_write_singularity()

    if os.path.isfile(TEMP_DOCKERFILE):
        os.remove(TEMP_DOCKERFILE)