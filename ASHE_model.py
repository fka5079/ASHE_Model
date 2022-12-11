# -*- coding: utf-8 -*-

# Import necessary libraries
import random
from random import randint
import numpy as np
import pandas as pd
import math
from scipy import interpolate
import re

class Workers:

    def __init__(self, time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process):
        # Each multiplier will be in the form of a list
        self.time = time
        self.takt_t = takt_t
        self.stress = stress
        self.complexity = complexity
        self.experience = experience
        self.procedures = procedures
        self.ergonomics = ergonomics
        self.FOD = FOD
        self.process = process
        self.time_mult = []
        self.complexity_mult = []
        self.hep = []
        
        # Interpolation to calculate time multiplier
        interp_col1 = np.asarray([self.takt_t, 2*self.takt_t, 5*self.takt_t, 50*self.takt_t])
        interp_col2 = np.asarray([10, 1, 0.1, 0.01])
        
        tck = interpolate.splrep(interp_col1, interp_col2, s=0, k=3)
        
        # Interpolation to calculate complexity multiplier
        interp_col3 = (np.asarray([0, 0.67, 1.33, 2.0]))
        interp_col4 = np.asarray([0.1, 1.0, 2.0, 5.0])
        
        for step in range(0, len(DM_mat)-1):
            self.time_mult.append(interpolate.splev(self.time[step], tck, der=0))
            self.complexity_mult.append(np.interp(self.complexity[step], interp_col3, interp_col4))
            psf = self.time_mult[step] * self.stress[step] * self.complexity[step] * self.experience[step] * self.procedures[step] * self.ergonomics[step] * self.FOD[step] * self.process[step]
            self.hep.append((0.01 * psf) / (0.01 * (psf - 1) + 1))


class Task:
    
    attempt = 0
    errorcount = 0
    
    # Status Key
    stat_complete = 4
    # stat_error_falsepos = 3
    stat_error_dependent = 2
    stat_error = 1
    stat_default = 0
    
    def __init__(self, DSM, bookshelf, errprob, a, b, c):
        self.DSM = DSM
        self.bookshelf = bookshelf
        self.errprob = errprob
        self.a = a  # deltat calculation constants
        self.b = b  # deltat calculation constants
        self.c = c  # deltat calculation constants
        self.steps = len(self.DSM) - 1
        
        # Some variables and lists to keep track
        error = [0, 1]
        falsepos = [0, 1]
        falsepos_prob = 0.01  # Probability of a false positive error occuring
        
            
        for step in range(1, self.steps + 1):
                
            deltat_tot = []  # A list that helps sum deltat at the end of a step
            # Resetting counters to track attempts and errors
            Task.attempt = 0
            Task.errorcount = 0
            
            # A Function that simulates a worker attempting a step.
            def Attempt(partx):
                
                self.part = partx
                # HEP when a step is being reattempted after an error is detected
                recurrant_errprob = self.errprob[step - 1]/2
                
                # Create a list of pre and post dependent parts
                depend = []  # List of all connected parts
                for y in range(1, len(self.DSM)):
                    if self.DSM[step, y] == "1":
                        depend.append(self.DSM[0, y])
                for x, dep_part in enumerate(depend):
                    if dep_part == self.part:
                        # Identify connected parts prior to current step:
                        self.depend_pre = depend[0 : x]
                        # Identify connected parts after current step:        
                        self.depend_post = depend[x+1 : len(depend)]
                        
                # Resetting counters to track attempts and errors
                #Task.attempt = 0
                #Task.errorcount = 0
                
                print(f"Attempting Step{step}")
                if self.bookshelf[self.part] == Task.stat_default:
                    
                    Task.attempt += 1
                    
                    # Error status in current step based on HEP for part
                    self.error = random.choices( error, weights = ((1-self.errprob[step - 1]), self.errprob[step - 1]), k=1 )
                    
                    # Probability of detecting error
                    self.errdetect = randint(0, 1)
                    
                    # If no error occurs -> mark part as complete
                    if self.error[0] == 0:
                        self.bookshelf[self.part] = Task.stat_complete
                        print(f"Step{step} was completed with no errors in " + str(Task.attempt) + " attempt(s).")
                        
                    # If an error occurs but is not detected -> update status of current and all post connected steps
                    # Use "try" because the last step in a task will not have any post dependent steps
                    elif self.error[0] == 1 and self.errdetect == 0:
                        Task.errorcount += 1
                        # Update deltat due to error that occured during first attempt
                        self.errorintensity = randint(1, 100)
                        # Calculate change in time based on error intensity
                        deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                        deltat_tot.append(deltat)
                        
                        try:
                            for x, dep_part in enumerate(self.depend_post):
                                self.bookshelf[dep_part] = Task.stat_error_dependent  # Update staus of post connected steps
                            print(f"Step{step} completed with un-detected error in " + str(Task.attempt) + " attempt(s).")
                        except:
                            # Last step has no post dependent steps -> add to error count and update current step only
                            print(f"Step{step} completed with un-detected error in " + str(Task.attempt) + " attempt(s).")
                        self.bookshelf[self.part] = Task.stat_error
                        
                    # If error occurs and is detected -> try to solve until no error (add to attempt and error count each time)
                    # Update time multiplier and error probability based on intensity of error
                    elif self.error[0] == 1 and self.errdetect == 1:
                        Task.errorcount += 1
                        # Update deltat due to error that occured during first attempt
                        self.errorintensity = randint(1, 100)
                        # Calculate change in time based on error intensity
                        deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                        deltat_tot.append(deltat)
                            
                        # Possibility that the error that occured is a false positive
                        self.falsepositive = random.choices( falsepos, weights = ((1-falsepos_prob), falsepos_prob), k=1 )
                        if self.falsepositive == 1:
                            print("A falsepositive error occured. Part was reattempted.")
                            Task.attempt += 1
                            Task.errorcount += 1
                            
                        elif self.falsepositive == 0:
                            # Calculate the intensity of the error
                            self.errorintensity = randint(1, 100)
                            # Calculate change in time based on error intensity
                            deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                            deltat_tot.append(deltat)
                            print(f"Error was detected in Step{step}.")
                            # Re-attempt task until detected error is solved
                            while self.error[0] == 1:
                                self.error = random.choices( error, weights = (( 1 - recurrant_errprob ), recurrant_errprob), k=1 )
                                # Calculate the intensity of the new error
                                self.errorintensity = randint(1, 100)
                                # Calculate change in time based on error intensity
                                deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                                deltat_tot.append(deltat)
                                Task.attempt += 1
                                Task.errorcount += 1
                                
                        self.bookshelf[self.part] = Task.stat_complete
                        print(f"Step{step} was completed in " + str(Task.attempt) + " attempt(s).")
           
                elif self.bookshelf[self.part] == Task.stat_error_dependent:
                    # Increase error counter because first attempt at step revealed dependency issue
                    Task.errorcount += 1
                    Task.attempt += 1
                    print(f"Dependency error detcted in Step{step}. Checking previous connected parts.")
                    # Calculate the intensity of the error
                    self.errorintensity = randint(1, 100)
                    # Calculate change in time based on error intensity
                    deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                    deltat_tot.append(deltat)
                    
                    # A function that is called when a part has a dependency error.
                    def Dependency(or_part, depend_prelist, depend_postlist):
                        or_step = int(re.findall("\d+", str(or_part))[0])
                        
                        stat_track = []
                        for x, part_1r in enumerate(depend_prelist):
                            step_1r = int(re.findall("\d+", str(part_1r))[0])  # Gives the number for current pre-dependent step & pep = prev erraneous part
                    
                            depend_1r = []
                            # List of parts dependent on part_1r
                            for y in range(1, len(self.DSM)):
                                if self.DSM[step_1r, y] == "1":
                                    depend_1r.append(self.DSM[0, y])
                            for x, dep_part in enumerate(depend_1r):
                                if dep_part == part_1r:
                                    # Identify connected parts prior to current step:
                                    depend_pre_1r = depend_1r[0 : x]
                                    # Identify connected parts after current step:        
                                    depend_post_1r = depend_1r[x+1 : len(depend_1r)]
                                    print(depend_post_1r)
                            
                            # When part_1r has a dependency error
                            if self.bookshelf[part_1r] == Task.stat_error_dependent:
                                # Attempt Part_1r but it cannot be solved until root error is solved.
                                try:
                                    attempts[step_1r - 1] += 1
                                    error_list[step_1r - 1] += 1
                                except:
                                    pass
                                
                                # Calculate deltat due to error in part_1r.
                                self.errorintensity = randint(1, 100)
                                # Calculate change in time based on error intensity
                                deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                                deltat_tot.append(deltat)
                                
                                print(f"{part_1r} has a dependency error. Checking previous connected parts.")
                                Dependency(part_1r, depend_pre_1r, depend_post_1r)
                                
                            # When part_1r has an error
                            elif self.bookshelf[part_1r] == Task.stat_error:
                                                                
                                # There is a chance for whether the error in part_1r is detected or not.
                                self.errdetect_1r = randint(0, 1)
                            
                                # When the error in part_1r is not detected, the worker moves on to the next dependent part.
                                if self.errdetect_1r == 0:
                                    print(f"{part_1r} has an error that was not detected. Attempting next part.")
                                
                                # When the error in part_1r is detected
                                elif self.errdetect_1r == 1:
                                    print(f"Error found in {part_1r}. Re-attempting {part_1r}.")
                                    attempts[step_1r - 1] += 1
                                    
                                    # Update hep of step_1r based on time available for step minus current time elapsed.
                                    time[step_1r - 1] = time[or_step - 1] - sum(deltat_tot)
                                    Worker = Workers(time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process)
                                    self.errprob = Worker.hep
                                    
                                    self.error_1r = random.choices( error, weights = ((1 - self.errprob[step_1r - 1]), self.errprob[step_1r - 1]), k=1 )
                                    
                                    # If the previous error in part_1r is solved
                                    if self.error_1r[0] == 0:
                                        print(f"Error in {part_1r} was solved.")
                                        self.bookshelf[part_1r] == Task.stat_complete
                                    
                                    # If the previous error in part_1r is not solved
                                    elif self.error_1r[0] == 1:
                                        print(f"Error reoccured in {part_1r}. Re-attempting {part_1r}.")
                                        error_list[step_1r - 1] += 1
                                        
                                        # Calculate deltat due to error in part_1r.
                                        self.errorintensity = randint(1, 100)
                                        # Calculate change in time based on error intensity
                                        deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                                        deltat_tot.append(deltat)
                                        
                                        # Update hep of step_1r based on time available for step minus current time elapsed.
                                        time[step_1r - 1] = time[or_step - 1] - sum(deltat_tot)
                                        Worker = Workers(time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process)
                                        self.errprob = Worker.hep
                                    
                                        self.error_1r = random.choices( error, weights = ((1-self.errprob[step_1r - 1]), self.errprob[step_1r - 1]), k=1 )
                                        
                                        while self.error_1r[0] == 1:
                                            error_list[step_1r - 1] += 1
                                            attempts[step_1r - 1] += 1
                                            
                                            print(f"Error in Step{step_1r} re-occured. Re-attempting step.")
                                            
                                            # Update hep of step_1r based on time available for step minus current time elapsed.
                                            time[step_1r - 1] = time[or_step - 1] - sum(deltat_tot)
                                            Worker = Workers(time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process)
                                            self.errprob = Worker.hep
                                    
                                            self.error_1r = random.choices( error, weights = ((1-self.errprob[step_1r - 1]), self.errprob[step_1r - 1]), k=1 )
                                            
                                            # Calculate deltat due to error in part_1r.
                                            self.errorintensity = randint(1, 100)
                                            # Calculate change in time based on error intensity
                                            deltat = self.a * math.exp((-self.errorintensity/self.b) ** self.c)
                                            deltat_tot.append(deltat)
                                            
                                    self.bookshelf[part_1r] = Task.stat_complete

                                    print(f"Re-attempting parts that are post dependent on {part_1r}.")
                                    for dep_part in depend_post_1r:
                                        if dep_part == or_part:
                                            #print(depend_postlist)
                                            #print(dep_part)
                                            break
                                        elif dep_part != or_part:
                                            #print(depend_postlist)
                                            #print(dep_part)
                                            Attempt(dep_part)
                                    self.bookshelf[or_part] = Task.stat_default
                                    # Attempt(or_part)
                                    
                            elif self.bookshelf[part_1r] == Task.stat_complete:
                                pass
                            
                            # if any of the pre-dependent parts have an error the current part cannot be solved
                            stat_track.append(int(self.bookshelf[part_1r]))
                        
                        for stat in stat_track:
                            if stat == Task.stat_error or stat == Task.stat_error_dependent:
                                print("Either an error or a dependency error in previous connected parts was not solved.")
                                self.bookshelf[or_part] = Task.stat_error_dependent
                                for dep_part in depend_postlist:
                                    self.bookshelf[dep_part] = Task.stat_error_dependent
                                break
                            elif stat != Task.stat_error and stat != Task.stat_error_dependent:
                                self.bookshelf[or_part] = Task.stat_complete
                        
                        # Now that previous connected parts have been checked, re-attempt current part
                        # This might take the worker back into dependency error checking
                        # Or result in the original part being solved
                        Attempt(or_part)

                    Dependency(self.part, self.depend_pre, self.depend_post)
                
                
            Attempt(partx = f"part{step}")
                    
            attempts.append(int(Task.attempt))
            error_list.append(int(Task.errorcount))
            print(attempts)
             
            # Add all delta_t values saved for this step
            # This constant gives the total change in time that will affect the next step and can be used to calculate new available time
            deltat_sum = sum(deltat_tot)
            # Update the time input for multiplier class; this works for all steps except the last step
            if step < (self.steps):
                time[step] -= deltat_sum    # Subtract time lost in error from next step
            elif step == (self.steps):
                pass
            # Use the new time list and call the multiplier class again to update HEP list
            Worker = Workers(time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process)
            self.errprob = Worker.hep


"""
Initiating class
"""

DSM = pd.read_csv("ball shooter_DSM.csv", header = None)
DM_mat = np.matrix(DSM)

# Replace header row and column with part numbers
parts = ['']
for part in range(1, len(DM_mat)):
    parts.append(f"part{part}")
DM_mat[0] = parts
DM_mat[:,0] = np.asarray([parts]).T

# Create a second version of the DSM with 0s as diagonals and no headers
DM_mat_mod = DM_mat

DM_mat_mod = np.delete(DM_mat_mod, 0, 0)
DM_mat_mod = np.delete(DM_mat_mod, 0, 1)

for N in range(0, len(DM_mat_mod)):
    DM_mat_mod[N, N] = 0
    
DM_mat_mod = DM_mat_mod.astype(int)

# Call Multipliers class and define values for each multiplier for each step within the task
# Input by user/shift supervisor
time = []  # Available time for each task
takt_t = 25  # Time needed for each task
stress = []
complexity = []
experience = []
procedures = []
ergonomics = []
FOD = []
process = []

# Create empty matrices that can be appended for results
attempts_mat = np.empty((0, len(DM_mat)-1), int)
# reattempts_mat = np.empty((0, len(DM_mat)-1), int)
error_list_mat = np.empty((0, len(DM_mat)-1), int)

# Constants for delta_t calculation
a = 1
b = 80     # scale parameter
c = 2      # shape parameter

# Run DSM through model x number of times
k = 0
while k < 10000:
    # Create dictionary to track errors
    bookshelf = { "part1":0, "part2":0, "part3":0, "part4":0, "part5":0, 
                  "part6":0, "part7":0, "part8":0, "part9":0, "part10":0, 
                  "part11":0, "part12":0, "part13":0, "part14":0, "part15":0, 
                  "part16":0, "part17":0, "part18":0, "part19":0, "part20":0,
                  "part21":0, "part22":0, "part23":0,  "part24":0, "part25":0, "part26":0,
                  "part27":0, "part28":0, "part29":0, "part30":0, "part31":0,
                  "part32":0, "part33":0, "part34":0, "part35":0, "part36":0, "part37":0,
                  "part38":0, "part39":0, "part40":0, "part41":0, "part42":0, "part43":0,
                  "part44":0, "part45":0, "part46":0, "part47":0, "part48":0, "part49":0,
                  "part50":0, "part51":0, "part52":0, "part53":0, "part54":0, "part55":0,
                  "part56":0, "part57":0, "part58":0, "part59":0, "part60":0, "part61":0,
                  "part62":0, "part63":0, "part64":0, "part65":0, "part66":0, "part67":0,
                  "part68":0, "part69":0, "part70":0, "part71":0, "part72":0, "part73":0,
                  "part74":0, "part75":0}
    
    #for item in parts[1:len(parts)]:
    #    bookshelf[item] = 0
    
    # Call Multipliers class and define values for each multiplier for each step within the task
    # Input by user/shift supervisor
    time = []  # Available time for each task
    takt_t = 25  # Time needed for each task
    stress = []
    complexity = []
    experience = []
    procedures = []
    ergonomics = []
    FOD = []
    process = []
    
    # Use following for loop to automatically fill a list with identical multiplier values
    # This section of the code can be updated to include a list with different values
    for item in range(0, len(DM_mat)-1):
        time.append(35)
        stress.append(1)
        # complexity.append(2)
        experience.append(0.5)
        procedures.append(1)
        ergonomics.append(1)
        FOD.append(1)
        process.append(1)
        
    # Calculating complexity for each step using in and out degree
    for ste_no in range(0, len(DM_mat_mod)):
        deg1 = np.sum(DM_mat_mod[ste_no,:])
        deg2 = np.sum(DM_mat_mod[:, ste_no])
        complexity.append(2.0*((deg1+deg2)/len(DM_mat_mod)))
    
    Worker = Workers(time, takt_t, stress, complexity, experience, procedures, ergonomics, FOD, process)
    # Errprob should be dependent on multipliers
    errprob = Worker.hep
    
    attempts = []
    error_list = []
    #reattempts = []
    #m = 0
    #while m < len(DM_mat) - 1:
        #reattempts.append(0)
        #m += 1
    
    Task1 = Task(DM_mat, bookshelf, errprob, a, b, c)
    
    attempts_mat = np.append(attempts_mat, np.array([attempts]), axis=0)
    # reattempts_mat = np.append(reattempts_mat, np.array([reattempts]), axis=0)
    error_list_mat = np.append(error_list_mat, np.array([error_list]), axis=0)
    print(bookshelf)
    print(f"Repetition {k+1} completed.")
    k += 1

# Save number of Attempts and Errors for each step, within each repetition
attempts_df = pd.DataFrame(attempts_mat)
# reattempts_df = pd.DataFrame(reattempts_mat)
error_list_df = pd.DataFrame(error_list_mat)
attempts_df.to_csv('Number_of_Attempts.csv')
# reattempts_df.to_csv('Number_of_ReAttempts.csv')
error_list_df.to_csv('Number_of_Errors.csv')