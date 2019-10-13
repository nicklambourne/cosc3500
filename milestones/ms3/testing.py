import subprocess
import multiprocessing as mp
import pandas as pd
from argparse import ArgumentParser
from os import mkdir, path
from random import choice
from sys import argv
from string import ascii_letters, digits, punctuation
from typing import List
from datetime import datetime
from matplotlib import pyplot
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from itertools import product
from shutil import rmtree


BINARY_PATH = "./bin/lcs"
AUTO_TEST_PATH = "./auto_test"
TEST_ROOT = "./test/"
SLURM_TEMPLATE = "./template.sh"

MAX_NODES = 4
MAX_TASKS = 64
MAX_TASKS_PER_NODE = 16
MAX_CPUS_PER_TASK = 64
MAX_CORES = 64


def write_template(output_path: str, template_file: str, context: dict) -> None:
    j2env = Environment(loader=FileSystemLoader(TEST_ROOT))
    template = j2env.get_template(template_file)
    content = template.render(context)
    with open(output_path, "w") as file_:
        file.write(content)


def b2_options(limit: int) -> List:
    return list([x ** 2 for x in range(1, limit + 1)])


class Job():
    def __init__(nodes: int, 
                 ntasks: int,
                 ntasks_per_node: int,
                 cpus_per_task: int):
        self.name = f"{name}-{nnodes}-{ntasks}-{ntasks_per_node}-{cpus_per_task}"
        self.nnodes = nodes
        self.ntasks = ntasks
        self.ntasks_per_node = ntasks_per_node
        self.cpus_per_task = cpus_per_task

    def generate_script():
        pass

    def run(n: int):
        subprocess.run("sh")


def produce_jobs() -> List[Job]:
    nodes = b2_options(MAX_NODES)
    ntasks = b2_options(MAX_TASKS)
    ntasks_per_node = b2_options(MAX_TASKS_PER_NODE)
    cpus_per_task = b2_options(MAX_CPUS_PER_TASK)
    
    combinations = product(nodes, ntasks, ntasks_per_node, cpus_per_task)
    valid_combinations = filter(lambda x: x[0] == x[1] * x[2] and x[1] * x[2] * x[3] <= MAX_CORES,
                                combinations)
    
    jobs = list()
    for combination in valid_combinations:
        job = Job(combination[0],
                  combination[1],
                  combination[2],
                  combination[3])
        jobs.append(job)

    return jobs


def random_string(length: int) -> str:
    return "".join(choice(ascii_letters+digits+punctuation) for i in range(length))


def run_single_test() -> None:
    pass


def generate_test_files(num_tests: int,
                        max_length: int) -> None:
    if path.isdir(TEST_PATH):
        rmtree(TEST_PATH)
    mkdir(TEST_PATH)
    for i in range(1, num_tests + 1):
        with open(f"{TEST_PATH}/test-{i}.txt", "w") as test_file:
            test_file.write(f"{random_string(int(i*(max_length/num_tests)))}\n"
                    f"{random_string(int(i*(max_length/num_tests)))}")


if __name__ == "__main__":
    jobs = produce_jobs()
    print(len(jobs))