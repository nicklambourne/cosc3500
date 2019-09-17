import multiprocessing as mp

from argparse import ArgumentParser
from os import mkdir
from random import choice
from string import printable
from typing import List


BINARY_PATH = "./bin/ass1"
TEST_PATH = "./test"


def random_string(length: int) -> str:
    return "".join(choice(printable) for i in range(length))


def run_single_test() -> None:
    pass

def generate_random_strings(num_strings: int) -> List[str]
    pool = Pool(mp.cpu_count() - 1)
    with pool as p:
        x = p.map(random_string)

def generate_test_files(num_tests: int,
                        max_length: int) -> None:
    os.mkdir(TEST_PATH)



if __name__ == "__main__":
    if len(argv) != 3:
        print("Usage: python test.py <num_tests> <max_length>")
        exit()
    parser = ArgumentParser()
    parser.add_argument("num_tests", type=int, help="Number of tests to run")
    parser.add_argument("max_length", type=int, help="Maximum length of strings to compare")
    generate_test_files(num_tests, max_length)
