import subprocess
import multiprocessing as mp
import pandas as pd
from argparse import ArgumentParser, Namespace
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
from math import ceil, log


BINARY_PATH = path.abspath("./bin/lcs-hybrid")
INPUT_ROOT = path.abspath("./test")
TEST_ROOT = path.abspath("./tests/")
TEMPLATE_ROOT = getcwd()
SLURM_TEMPLATE = "template.sh"
INPUT_NAME = "xlong_in.txt"

MAX_NODES = 4
MAX_TASKS = 64
MAX_TASKS_PER_NODE = 16
MAX_CPUS_PER_TASK = 16
MAX_CORES_PER_NODE = 16


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
                 cpus_per_task: int,
                 input_file: str = path.abspath(path.join(INPUT_ROOT, INPUT_NAME)),
                 output_file: str = None):
        self.name = f"{nnodes}-{ntasks}-{ntasks_per_node}-{cpus_per_task}"
        self.directory = path.abspath(path.join(TEST_ROOT, self.name))
        self.script = path.abspath(path.join(self.directory, "go.sh"))
        self.input = input_file
        self.output = output_file if output_file else path.abspath(path.join(self.directory, self.name + ".out"))
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
        return f"{self.nnodes}-{self.ntasks}-{self.ntasks_per_node}-{self.cpus_per_task}-{self.input}"

def produce_strong_jobs() -> List[Job]:
    nnodes = b2_options(MAX_NODES)
    ntasks = b2_options(MAX_TASKS)
    ntasks_per_node = b2_options(MAX_TASKS_PER_NODE)
    cpus_per_task = b2_options(MAX_CPUS_PER_TASK)
    
    combinations = product(nnodes, ntasks, ntasks_per_node, cpus_per_task)
    valid_combinations = filter(lambda x: x[1] == x[0] * x[2] and 
                                (x[1] * x[3] <= MAX_CORES_PER_NODE or x[3] == 1),
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
    return "".join(choice(ascii_letters+digits) for i in range(length))


def write_weak_file(test_string_a: str, test_string_b: str, length: int):
    with open(f"{INPUT_ROOT}/{length}.in", "w") as input_file:
        input_file.write(f"{test_string_a[:length]}\n"
                            f"{test_string_b[:length]}")



def produce_weak_jobs() -> List[Job]:
    increments = [10, 100, 1000]
    test_string_a = random_string(MAX_TASKS * max(increments))
    test_string_b = random_string(MAX_TASKS * max(increments))
    
    jobs = list()

    for increment in increments:
        for ntasks in b2_options(MAX_TASKS):
            length = increment * ntasks
            write_weak_file(test_string_a, test_string_b, length)
            nnodes = ceil(ntasks / MAX_TASKS_PER_NODE)
            ntasks_per_node = ntasks if ntasks <= MAX_TASKS_PER_NODE else int(ntasks / nnodes)
            cpus_per_task = 1
            job = Job(nnodes, ntasks, ntasks_per_node, cpus_per_task, 
                      f"{INPUT_ROOT}/{length}.in", f"{INPUT_ROOT}/{length}.out")
            jobs.append(job)
    
    return jobs


def generate_test_files(num_tests: int,
                        max_length: int) -> None:
    if path.isdir(TEST_PATH):
        rmtree(TEST_PATH)
    mkdir(TEST_PATH)
    for i in range(1, num_tests + 1):
        with open(f"{TEST_PATH}/test-{i}.txt", "w") as test_file:
            test_file.write(f"{random_string(int(i*(max_length/num_tests)))}\n"
                    f"{random_string(int(i*(max_length/num_tests)))}")


def run_jobs(args: Namespace) -> None:
    jobs = None
    if args.test_type == "strong":
        jobs = produce_strong_jobs()
    elif args.test_type == "weak":
        jobs = produce_weak_jobs()
    
    for job in jobs:
        print(f"{job.name}: ", end="", flush=True)
        if not args.dry:
            job.run(args.num_tests)
        else:
            print(flush=True)

    if args.dry:
        print(f"\nEnd dry run of {len(jobs)} jobs")
    else:
        print(f"\nEnd run of {len(jobs)} jobs")


def run_report(args: Namespace) -> None:
    pass


if __name__ == "__main__":
    if not path.exists(TEST_ROOT):
        mkdir(TEST_ROOT)

    master_parser = ArgumentParser(usage="python -m testing <subcommand> [options]")

    subcommand_parser = master_parser.add_subparsers(dest="command")
    run_parser = subcommand_parser.add_parser("run", help="Run tests using slurm")
    run_parser.add_argument("--dry", dest="dry", action="store_true", 
                            help="A dry run will not launch jobs")
    run_parser.add_argument("test_type", choices=["strong", "weak"],
                            help="Whether to run tests for strong or weak scaling")
    run_parser.add_argument("num_tests", type=int, help="Number of tests to run")
    
    report_parser = subcommand_parser.add_parser("report", help="Produce report from test results")
    report_parser.add_argument("test_folder", help="Directory to run reporting on")

    args = master_parser.parse_args()

    if args.command == "run":
        run_jobs(args)
    elif args.command == "report":
        run_report(args)
    else:
        pass
