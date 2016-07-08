Running JED is now greatly simplified!

You run the program using executable jar files:

JED_Driver.jar    and   JED_Batch_Driver.jar   are the two executable java jar files.

The Driver is the 'singleton' variant while the Batch_Driver allows you to specify multiple jobs.

Simply copy the jar files to your home directory, for example, put them in ~/java/bin

TO RUN THE PROGRAM USE THE FOLLOWING COMMANDS:

java -d64 -jar JED_Driver.jar
where it is assumed that your input file is named JED_Driver.txt and is in the same directory as JED_Driver.jar

java -d64 -jar JED_Driver.jar "/path/to/your/input/file/JED_Driver.txt"
where you supply the path to your input file as an argument to the command.
Note: You can name your input file as you like when doing this ;-)

**********************************************************************************************************************************************************

java -d64 -jar JED_Batch_Driver.jar
where it is assumed that your input file is named 'JED_Batch_Driver.txt' and is in the same directory as JED_Batch_Driver.jar


java -d64 -jar JED_Batch_Driver.jar "/path/to/your/input/file/JED_Batch_Driver.txt"
where you supply the path to your input file as an argument to the command.
Note: You can name your input file as you like when doing this ;-)


SAMPLE PBS SCRIPTS ARE IN THE /dcm/JED directory, for single jobs and using job arrays.

NOTE: THE BATCH DRIVER IS MEANT FOR RUNNING JED JOBS SERIALLY, AS WHEN YOU DO NOT HAVE CLUSTER RESOURCES AVAILABLE.
          IF YOU HAVE CLUSTER ACCESS, THEN USE JED DRIVER AND A PBS SCRIPT THAT IMPLEMENTS A JOB ARRAY.

THE POOL DATA AND SUBSPACE ANALYSIS DRIVERS WORK THE SAME WAY.

Have Fun!!!