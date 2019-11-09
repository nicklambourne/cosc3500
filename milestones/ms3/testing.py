import subprocess
import multiprocessing
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser, Namespace
from os import mkdir, path, chdir, getcwd
from random import choice
from sys import argv
from string import ascii_letters, digits, punctuation
from typing import List
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from itertools import product
from shutil import rmtree
from math import ceil, log
from os import walk
from uuid import uuid4


BINARY_PATH = path.abspath("./bin/lcs-hybrid")
INPUT_ROOT = path.abspath("./test")
TEST_ROOT = path.abspath("./tests/")
TEMPLATE_ROOT = getcwd()
SLURM_TEMPLATE = "template.sh"
INPUT_FILE_NAME = path.abspath(f"./{str(uuid4())}.txt")

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
                 input_file: str = path.abspath(path.join(INPUT_ROOT, INPUT_FILE_NAME)),
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


def produce_strong_jobs(args: Namespace) -> List[Job]:
    with open(INPUT_FILE_NAME, "w") as input_file:
        input_file.write(f"{random_string(args.input_size)}\n"
                         f"{random_string(args.input_size)}\n")
    
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
                  combination[3],
                  f"{INPUT_FILE_NAME}", 
                  f"{INPUT_ROOT}/{str(uuid4())}.out")
        jobs.append(job)

    return jobs


def random_string(length: int) -> str:
    return "".join(choice(ascii_letters+digits) for i in range(length))


def write_weak_file(test_string_a: str, test_string_b: str, length: int) -> None:
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
            nnodes = ceil(ntasks / MAX_TASKS_PER_NODE)
            ntasks_per_node = ntasks if ntasks <= MAX_TASKS_PER_NODE else int(ntasks / nnodes)
            for cpus_per_task in [1, 2, 4, 8, 16]: 
                if ntasks_per_node * cpus_per_task <= MAX_CORES_PER_NODE:
                    length = increment * ntasks * cpus_per_task
                    write_weak_file(test_string_a, test_string_b, length)
                    job = Job(nnodes, ntasks, ntasks_per_node, cpus_per_task, 
                              f"{INPUT_ROOT}/{length}.in", f"{INPUT_ROOT}/{length}.out")
                    jobs.append(job)
    
    return jobs


def run_jobs(args: Namespace) -> None:
    jobs = None
    if args.test_type == "strong":
        jobs = produce_strong_jobs(args)
    elif args.test_type == "weak":
        jobs = produce_weak_jobs()
    
    for job in jobs:
        print(f"{job.name}: ", end="", flush=True)
        if not args.dry:
            job.run(args.num_tests)
        else:
            print(flush=True)

    print(f"\nEnd{' dry' if args.dry else ''} run of {len(jobs)} jobs and {len(jobs) * args.num_tests} passes\n")


def run_report(args: Namespace) -> None:
    results = list()
    for root, dirs, files in walk(args.test_folder):
        for file in files:
            if file.endswith('.out'):
                results.append(f"{root}/{file}")

    findings = list()

    for result in results:
        with open(result, "r") as result_file:
            try:  # this allows us to run over a running directory and skip malformed files
                settings_line = result_file.readline()
                setting_components = settings_line.split(" ")[0].split("-")
                nnodes = int(setting_components[0])
                ntasks = int(setting_components[1])
                ntasks_per_node = int(setting_components[2])
                cpus_per_task = int(setting_components[3])

                size_line = result_file.readline()
                problem_size = int(size_line[1:-2].split(",")[0])

                result_file.readline()  # Throw this away -> output string

                run_time = float(result_file.readline().split(" ")[0])

                findings.append([nnodes, ntasks, ntasks_per_node, cpus_per_task,
                                 problem_size, int(problem_size/(ntasks)), run_time])

            except:
                pass
    
    dataframe = pd.DataFrame(columns=["nnodes", "ntasks", "ntasks_per_node", 
                                      "cpus_per_task", "problem_size", "magnitude",
                                      "run_time"],
                            data=findings)
        
    dataframe = dataframe.sort_values(["nnodes", "ntasks", "cpus_per_task"])

    if args.test_type == "strong": 
        report_dataframe = dataframe.groupby(["nnodes", "ntasks", "ntasks_per_node", "cpus_per_task"]).mean()
        report_dataframe.to_csv(f"{args.test_type}.csv")
        dataframe = dataframe.groupby(["ntasks", "cpus_per_task"]).mean()
        print(report_dataframe)
        
        
    elif args.test_type == "weak": 
        dataframe = dataframe.groupby(["nnodes", "ntasks", "ntasks_per_node", 
                        "cpus_per_task", "problem_size", "magnitude"]).mean()
        dataframe = dataframe.sort_values(["magnitude", "ntasks"])

    dataframe.reset_index(inplace=True)

    figure, axis = plt.subplots()
    
    if args.test_type == "weak":
        for key, group in dataframe.groupby("magnitude"):
            magnitude = group.magnitude.unique()[0]
            axis = group.plot(kind="line", x="ntasks", y="run_time", label=key, linestyle="-", marker="o")
            plt.legend(title="Size/Processing Unit", loc="best")
            axis.set_ylabel("Wall Time (s)")
            axis.set_xlabel("Processing Units (ntasks * cpus_per_task)")
            plt.savefig(f"weak-{magnitude}.png")
        dataframe.to_csv(f"{args.test_type}.csv")
        print(dataframe)

    elif args.test_type == "strong":
        dataframe = dataframe.sort_values(["ntasks", "cpus_per_task"])
        index = 0
        colours = ["red", "blue", "orange", "green", "purple"]
        plt.xscale('log', basex=2)
        for key, group in dataframe.groupby("cpus_per_task"):
            axis = group.plot(ax=axis, kind="line", x="ntasks", y="run_time", label=key, c=colours[index], linestyle="-", marker="o")
            index += 1
        plt.legend(title="Thread Count", loc="best")
        axis.set_ylabel("Wall Time (s)")
        axis.set_xlabel("MPI Process Count (ntasks)")
        axis.set_xticklabels(["0", "1", "2", "4", "8", "16", "32", "64"])
        plt.savefig(f"strong.png")


if __name__ == "__main__":
    if not path.exists(TEST_ROOT):
        mkdir(TEST_ROOT)

    master_parser = ArgumentParser(usage="python -m testing <subcommand> [options]")

    subcommand_parser = master_parser.add_subparsers(dest="command")

    run_parser = subcommand_parser.add_parser("run", help="Run tests using slurm")
    run_parser.add_argument("--dry", dest="dry", action="store_true", 
                            help="A dry run will not launch jobs")

    run_type_parser = run_parser.add_subparsers(dest="test_type")

    strong_parser = run_type_parser.add_parser("strong",
                            help="Run tests for strong scaling")
    strong_parser.add_argument("input_size", type=int)
    strong_parser.add_argument("num_tests", type=int, help="Number of tests to run")
    weak_parser = run_type_parser.add_parser("weak",
                            help="Run tests for weak scaling")
    weak_parser.add_argument("num_tests", type=int, help="Number of tests to run")
    
    
    report_parser = subcommand_parser.add_parser("report", 
                                                 help="Produce report from test results")
    report_parser.add_argument("test_type", choices=["strong", "weak"],
                               help="Whether to run reporting for strong or weak scaling")
    report_parser.add_argument("test_folder", help="Directory to run reporting on")

    args = master_parser.parse_args()

    if args.command == "run":
        run_jobs(args)
    elif args.command == "report":
        run_report(args)
    else:
        pass
