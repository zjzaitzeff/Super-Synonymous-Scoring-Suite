This is my first documentation. I'm sorry...

Congradulations! You downloaded this off the github. Currently PITA, and Sei do not work. We hope to get PITA working soon and Sei will be switched out for Enformer.

DOCKER SETUP
------------------------------------
In order to run this pipeline you must have docker desktop installed and running on your computer.
You can find information on how to install docker here: https://www.docker.com/products/docker-desktop/
------------------------------------

PIPELINE SETUP
-------------------------------------
Once docker is running run in your terminal the following commands in order from the pipeline_docker_practice:
cd pipeline
docker build -t pipeline_snv .
        The step above should take around 3290.5s.

Once complete you can run the following examples in your terminal:

docker run --rm -v $(pwd)/../inputs_test:/pipeline/inputs_test -v $(pwd)/../more_testing:/pipeline/more_testing pipeline_snv /pipeline/inputs_test /pipeline/more_testing /pipeline/ref/miR_Family_Info_Human.fa 500

docker run --rm -v $(pwd)/../inputs_test:/pipeline/inputs_test -v $(pwd)/../outputs_test:/pipeline/outputs_test pipeline_snv /pipeline/inputs_test /pipeline/outputs_test /pipeline/ref/miR_Family_Info_Human.fa 500
------------------------------------

Important Notes
-------------------------------------
Inorder to generate output files onto your computer you must mount them unto the docker in the command line. Mounting should occur before "pipeline_snv" in the command line. The names "inputs_test" and "more_testing" are not hard coded into the pipeline so you can name them however you would like when you mount them. 

Ex: -v $(pwd)/../inputs_test:/pipeline/inputs_test or -v $(pwd)/../more_testing:/pipeline/more_testing


