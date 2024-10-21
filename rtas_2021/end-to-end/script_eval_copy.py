#!/usr/bin/env python3

import gc  # garbage collector
import argparse
import math
import numpy as np
import utilities.chain as c
import utilities.communication as comm
import utilities.generator_WATERS as waters
import utilities.generator_UUNIFAST as uunifast
import utilities.transformer as trans
import utilities.event_simulator as es
import utilities.analyzer as a
import utilities.evaluation as eva
import time


# Generate tasksets from WATERS benchmark
print("WATERS benchmark.")
r = 1000     # number of runs
u = 80.0    # utilization [%]


# Statistical distribution for task set generation from table 3
# of WATERS free benchmark paper.
profile = [0.03 / 0.85, 0.02 / 0.85, 0.02 / 0.85, 0.25 / 0.85,
            0.25 / 0.85, 0.03 / 0.85, 0.2 / 0.85, 0.01 / 0.85,
            0.04 / 0.85]
# Required utilization:
req_uti = 80.0/100.0
# Maximal difference between required utilization and actual
# utilization is set to 1 percent:
threshold = 1.0

# Create task sets from the generator.
# Each task is a dictionary.
print("\tCreate task sets.")
task_sets_waters = []
while len(task_sets_waters) < r:
    task_sets_gen = waters.gen_tasksets(
            1, req_uti, profile, True, threshold/100.0, 4)
    #print(task_sets_gen[0])
    task_sets_waters.append(task_sets_gen[0])

# Transform tasks to fit framework structure.
# Each task is an object of utilities.task.Task.
trans1 = trans.Transformer("1", task_sets_waters, 10000000)
task_sets = trans1.transform_tasks(False)

# Set an array for task counts for each taskset
ts_t_size = []
# Display and save generated tasksets
with open('taskset_80_1000_tmp.txt', 'w') as file:
    i = 0
    for task in task_sets:
        ts_t_size.append(len(task))
        print("Task set ", i, "(period, wcet, deadline, id, priority, phase, bcet)")
        #file.write(f"{'T':>10} {'C':>10} {'D':>10} {'ID':>5} {'Priority':>10} {'Phase':>10} {'BCET':>10}\n")
        #file.write(f"Task set {i} \n")
        for t in task:
            print(t.period, t.wcet, t.deadline, t.id, t.priority, t.phase, t.bcet)
            file.write(f"{t.period/10000000:>10} {t.wcet/10000000:>10} {t.deadline/10000000:>10} {t.id:>5} {t.priority:>10}\n")
            #print(task[0].id, task[0].phase, task[0].bcet, task[0].wcet)
        i += 1

# generate ce chains
ce_chains = waters.gen_ce_chains(task_sets)

# Set an array for chain counts for each taskset
ts_ch_size = []

# Display and save generated chains
with open('chains_80_1000_tmp.txt', 'w') as file:
    #file.write(f"number of chains : {len(ce_chains)}\n")
    for ce in ce_chains:    
        print("number of chains in a taskset :", len(ce))
        ts_ch_size.append(len(ce))
        #file.write(f"Length of ce : {len(ce)}\n")
        #file.write(f"chain  id  id  id  id  id  id\n")
        i = 1
        for c in ce:
            print("Chain ", i)
            file.write(f"{i:>5}")
            for task in c.chain:
                print("Task id is ", task.id)
                file.write(f"{task.id:>4}") 
            file.write(f"\n") 
            i += 1


###
# First analyses (TDA, Davare, Duerr).
###
debug_flag = True
print("=First analyses (TDA, Davare, Duerr).=")
analyzer = a.Analyzer("0")

# Save task-count and chain-count of each taskset 
with open('t_ch_each_ts_80_1000_tmp.txt', 'w') as file:
    file.write(f"{ts_t_size}\n")
    file.write(f"{ts_ch_size}\n")

try:
    # TDA for each task set.
    print("TDA.")
    for idxx in range(len(task_sets)):
        try:
            # TDA.
            i = 1
            for task in task_sets[idxx]:
                # Prevent WCET = 0 since the scheduler can
                # not handle this yet. This case can occur due to
                # rounding with the transformer.
                if task.wcet == 0:
                    raise ValueError("WCET == 0")
                task.rt = analyzer.tda(task, task_sets[idxx][:(i - 1)])
                if task.rt > task.deadline:
                    raise ValueError(
                            "TDA Result: WCRT bigger than deadline!")
                i += 1
        except ValueError:
            # If TDA fails, remove task and chain set and continue.
            task_sets.remove(task_sets[idxx])
            ce_chains.remove(ce_chains[idxx])
            continue

    # End-to-End Analyses.
    print("Test: Davare.")
    analyzer.davare(ce_chains)

    #print("Test: Duerr Reaction Time.")
    #analyzer.reaction_duerr(ce_chains)

    #print("Test: Duerr Data Age.")
    #analyzer.age_duerr(ce_chains)

    # Start the timer
    start_time = time.time()

    # Analysis RTAS 2021
    i = 0  # task set counter
    schedules = []
    # To save 'max data age', 'reduced-max data age', 'reaction time' as a file
    with open('reaction_t_80_1000_tmp.txt', 'w') as file: 
        #file.write(f"chn_id  react_t \n")
        for task_set in task_sets:
            print("=Task set ", i+1)
            #file.write(f"Task set : {i+1}\n") 

            # Skip if there is no corresponding cause-effect chain.
            if len(ce_chains[i]) == 0:
                continue

            # Event-based simulation.
            print("Simulation.")

            simulator = es.eventSimulator(task_set)

            # Determination of the variables used to compute the stop
            # condition of the simulation
            max_e2e_latency = max(ce_chains[i], key=lambda chain:
                                    chain.davare).davare
            max_phase = max(task_set, key=lambda task: task.phase).phase
            max_period = max(task_set, key=lambda task: task.period).period
            hyper_period = analyzer.determine_hyper_period(task_set)

            sched_interval = (
                    2 * hyper_period + max_phase  # interval from paper
                    + max_e2e_latency  # upper bound job chain length
                    + max_period)  # for convenience

            # Information for end user.
            print("\tNumber of tasks: ", len(task_set))
            print("\tHyperperiod: ", hyper_period)
            number_of_jobs = 0
            for task in task_set:
                number_of_jobs += sched_interval/task.period
            print("\tNumber of jobs to schedule: ",
                    "%.2f" % number_of_jobs)

            # Stop condition: Number of jobs of lowest priority task.
            simulator.dispatcher(
                    int(math.ceil(sched_interval/task_set[-1].period)))

            # Simulation without early completion.
            schedule = simulator.e2e_result()
            schedules.append(schedule)

            # Analyses.
            j = 1
            for chain in ce_chains[i]:
                #print("Test: Our Data Age.")
                #analyzer.max_age_our(schedule, task_set, chain, max_phase,
                #                        hyper_period, reduced=False)
                #analyzer.max_age_our(schedule, task_set, chain, max_phase,
                #                        hyper_period, reduced=True)

                print("Test: Our Reaction Time.")
                analyzer.reaction_our(schedule, task_set, chain, max_phase,
                                        hyper_period)

                #print("Max age : ", chain.our_age)
                #print("Max age (reduced) : ", chain.our_red_age)
                print("Max reaction : ", chain.our_react)

                file.write(f"{j:>5}")
                #file.write(f"{chain.our_age/10000000:>10.3f}")
                #file.write(f"{chain.our_red_age/10000000:>10.3f}")
                file.write(f"{chain.our_react/10000000:>10.3f}\n")
                j += 1

                # Kloda analysis, assuming synchronous releases.
                #print("Test: Kloda.")
                #analyzer.kloda(chain, hyper_period)            
            i += 1
    # End the timer
    end_time = time.time()
    # Calculate and print the elapsed time
    elapsed_time = end_time - start_time    
    
    # Open the file in append mode
    with open('t_ch_each_ts_80_1000_tmp.txt', 'a') as file:
        # Add the elapsed time at the end
        file.write(f"Elapsed time: {elapsed_time:.4f} seconds\n")
    print(f"Elapsed time: {elapsed_time:.4f} seconds")      

except Exception as e:
            print(e)
            print("ERROR: analysis")
            if debug_flag:
                breakpoint()
            else:
                task_sets = []
                ce_chains = []
