import subprocess
import multiprocessing as mp
import pandas as pd
from argparse import ArgumentParser
from os import mkdir, path, chdir, getcwd
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
from math import log


BINARY_PATH = path.abspath("./bin/lcs-hybrid")
AUTO_TEST_PATH = "./auto_test"
INPUT_ROOT = path.abspath("./test")
TEST_ROOT = path.abspath("./tests/")
TEMPLATE_ROOT = getcwd()
SLURM_TEMPLATE = "template.sh"
INPUT_NAME = "xlong_in.txt"

MAX_NODES = 4
MAX_TASKS = 64
MAX_TASKS_PER_NODE = 16
MAX_CPUS_PER_TASK = 16
MAX_CORES = 64


def write_template(output_path: str, template_file: str, context: dict) -> None:
    j2env = Environment(loader=FileSystemLoader(searchpath=TEMPLATE_ROOT))
    template = j2env.get_template(template_file)
    content = template.render(context)
    with open(output_path, "w") as file_:
        file_.write(content)


def b2_options(limit: int) -> List:
    return list([2 ** x for x in range(0, int(log(limit, 2) )+ 1)])


class Job():
    def __init__(self,
                 nnodes: int, 
                 ntasks: int,
                 ntasks_per_node: int,
                 cpus_per_task: int):
        self.name = f"{nnodes}-{ntasks}-{ntasks_per_node}-{cpus_per_task}"
        self.directory = path.abspath(path.join(TEST_ROOT, self.name))
        self.script = path.abspath(path.join(self.directory, "go.sh"))
        self.input = path.abspath(path.join(INPUT_ROOT, INPUT_NAME))
        self.output = path.abspath(path.join(self.directory, self.name + ".out"))
        self.binary = BINARY_PATH
        self.nnodes = nnodes
        self.ntasks = ntasks
        self.ntasks_per_node = ntasks_per_node
        self.cpus_per_task = cpus_per_task

    def run(self, n: int):
        if not path.exists(self.directory):
            mkdir(self.directory)

        context = {key: value for key, value in self.__dict__.items() \
                    if not key.startswith('__') and not callable(key)}

        write_template(self.script, SLURM_TEMPLATE, context)

        chdir(self.directory)

        for _ in range(n):
            subprocess.run(["sbatch", self.script])

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"{self.nnodes}-{self.ntasks}-{self.ntasks_per_node}-{self.cpus_per_task}"

def produce_jobs() -> List[Job]:
    nnodes = b2_options(MAX_NODES)
    ntasks = b2_options(MAX_TASKS)
    ntasks_per_node = b2_options(MAX_TASKS_PER_NODE)
    cpus_per_task = b2_options(MAX_CPUS_PER_TASK)
    
    combinations = product(nnodes, ntasks, ntasks_per_node, cpus_per_task)
    valid_combinations = filter(lambda x: x[1] == x[0] * x[2] and x[1] * x[3] <= MAX_CORES,
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
    if not path.exists(TEST_ROOT):
        mkdir(TEST_ROOT)
    
    if len(argv) != 2:
        print("Usage: python testing.py <num_tests>")
        exit()

    parser = ArgumentParser()
    parser.add_argument("num_tests", type=int, help="Number of tests to run")
    args = parser.parse_args()

    jobs = produce_jobs()
    
    for job in jobs:
        job.run(args.num_tests)