import multiprocessing as mp
from shutil import rmtree
import subprocess
import pandas as pd
from argparse import ArgumentParser
from os import mkdir, path
from random import choice
from sys import argv
from string import ascii_letters, digits, punctuation
from typing import List
from datetime import datetime
from matplotlib import pyplot


BINARY_PATH = "./bin/lcs-serial"
TEST_PATH = "./auto_test"


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

def run_tests(num_tests: int,
              max_length: int) -> pd.DataFrame:
    run_times = []
    for i in range(1, num_tests + 1):
        start_time = datetime.now()
        subprocess.run([BINARY_PATH, f"{TEST_PATH}/test-{i}.txt", f"{TEST_PATH}/test-{i}.out"])
        end_time = datetime.now()
        run_times.append((end_time - start_time).total_seconds())
    data = pd.DataFrame({"string_length": range(1, max_length + 1, int(max_length/num_tests)),
                         "run_times": run_times,
                         "test_no": range(1, num_tests + 1)})
    return data

def produce_graph(data: pd.DataFrame) -> None:
    scatter = data.plot.scatter(x="string_length", y="run_times")
    pyplot.xlabel("String Length (characters)")
    pyplot.ylabel("Run Time (wall - seconds)")
    pyplot.title("Run Time Performance of LCS")
    pyplot.grid(True)
    pyplot.savefig("runtime_performance.png", format="png")
    pyplot.show()

if __name__ == "__main__":
    if len(argv) != 3:
        print("Usage: python testing.py <num_tests> <max_length>")
        exit()
    parser = ArgumentParser()
    parser.add_argument("num_tests", type=int, help="Number of tests to run")
    parser.add_argument("max_length", type=int, help="Maximum length of strings to compare")
    args = parser.parse_args()
    generate_test_files(args.num_tests, args.max_length)
    data = run_tests(args.num_tests, args.max_length)
    produce_graph(data)

