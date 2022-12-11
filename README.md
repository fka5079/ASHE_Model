Hello! This documents provides a quick rundown on how to run the model included in this repository: An Agent-Based Model that Simulates Human Error (ASHE)

As the model stands currently, it includes two classes that represent a worker and a task.
These classes work with each other to utilize the inputs provided to the model and predict the liklihood of human error by running the simulation a given number of times.

INPUTS:
Design Structure Matrix:
A Design Structure Matrix (DSM) of choice is imported into the model at line 352.
A sample DSM, for the purposes of running this model has been included and is named "ball shooter_DSM.csv"
This matrix can be replaced by another DSM of any size.
If the DSM imported has a diagnoal of 1s, no changes need to be made to run the model.
However, if the DSM imported has a diagonal of 0s, lines 362 to 371 are not needed.

Constants:
Lines 391 to 393 can be used adjust constants a, b, and c, which dictate the calculation of change in time due to the occurence of an error (explained further in paper).

Model Repetitions:
Use line 397 to adjust the number of repetitions. Currently, the agent works through the given task 10,000 times.

SPAR-H Multipliers:
The next inputs needed to run ASHE include 7 Performance Shaping Factors (PSFs). Currently, PSFs are set to values of 1 or 0.5.
These values get updated as the model runs. However, their initial status can beadjusted using lines 419 yo 439

SAVING RESULTS:
Lines 469 to 474 can be used to adjust where the outputs of the model are saved.
Currently, the model outputs 2 csv files that track the number of attempts that the agent took to complete each step within the task for each repetition and the number of errors that occurred during the process.
The "ball shooter_DSM.csv" file was run through the ASHE model and sample outputs produced were included for reference.
Please keep in mind that this model is stochastic and the results will vary each time.

Thank you!
