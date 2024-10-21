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


# Generate tasksets from WATERS benchmark
print("WATERS benchmark.")
r = 1      # number of runs
u = 50.0    # utilization [%]


# Statistical distribution for task set generation from table 3
# of WATERS free benchmark paper.
profile = [0.03 / 0.85, 0.02 / 0.85, 0.02 / 0.85, 0.25 / 0.85,
            0.25 / 0.85, 0.03 / 0.85, 0.2 / 0.85, 0.01 / 0.85,
            0.04 / 0.85]
# Required utilization:
req_uti = 50.0/100.0
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

# Display generated tasksets
i = 0
for task in task_sets:
    print("Task set ", i, "(id, phase, bcet, wcet, period, deadline, priority)")
    for t in task:
        print(t.id, t.phase, t.bcet, t.wcet, t.period, t.deadline, t.priority)
        #print(task[0].id, task[0].phase, task[0].bcet, task[0].wcet)
    i += 1

# generate ce chains
ce_chains = waters.gen_ce_chains(task_sets)

# Display generated chains
print("number of chains : ", len(ce_chains))
for ce in ce_chains:    
    print("Length of ce :", len(ce))
    i = 0
    for c in ce:
        print("Chain ", i)
        for task in c.chain:
            print("Task id is ", task.id)
        i += 1




###
# First analyses (TDA, Davare, Duerr).
###
debug_flag = True
print("=First analyses (TDA, Davare, Duerr).=")
analyzer = a.Analyzer("0")

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

    # Analysis RTAS 2021
    i = 0  # task set counter
    schedules = []
    for task_set in task_sets:
        print("=Task set ", i+1)

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
        for chain in ce_chains[i]:
            print("Test: Our Data Age.")
            analyzer.max_age_our(schedule, task_set, chain, max_phase,
                                    hyper_period, reduced=False)
            analyzer.max_age_our(schedule, task_set, chain, max_phase,
                                    hyper_period, reduced=True)

            print("Test: Our Reaction Time.")
            analyzer.reaction_our(schedule, task_set, chain, max_phase,
                                    hyper_period)

            print("Max age : ", chain.our_age)
            print("Max age (reduced) : ", chain.our_red_age)
            print("Max reaction : ", chain.our_react)

            # Kloda analysis, assuming synchronous releases.
            #print("Test: Kloda.")
            #analyzer.kloda(chain, hyper_period)            
        i += 1

except Exception as e:
            print(e)
            print("ERROR: analysis")
            if debug_flag:
                breakpoint()
            else:
                task_sets = []
                ce_chains = []
