import argparse
import math
import numpy as np

# Two cases:
# - Case 1: Given starting numbers, print sequence
# - Case 2: Given limiting range of (a, a), enumerate instances of each cycle

# upper bound of the prime list we produce
MAX_PRIME = 1000000

# list needed for ordered factorization
prime_list = []

# dictionary of primes? (will implement bit array if needed)
prime_bools = []

# cache of nodes visited

class Node:

    def __init__(self):
        self.length = 0
        self.terms = []
        self.next = None
        self.cycle = None

    def index(self):
        return tuple(self.terms[:2])


nodes = {}


# class of cycles

class Cycle:

    def __init__(self):
        self.nodes = set()
        self.length = 0

cycles = {}

# cache of reductions
reductions = np.zeros(MAX_PRIME, dtype=np.uint64)

def generate_prime_lists(size):

    global prime_list, prime_bools

    prime_bools = np.ones(size, dtype=bool)
    prime_bools[0] = False
    prime_bools[1] = False

    # sieve of eratosthenes first
    index = 1

    # generate boolean by sieve of eratosthenes
    for k, v in enumerate(prime_bools):
        if not v:
            continue
        i = 2*k
        while i < len(prime_bools):
            prime_bools[i] = False 
            i += k

    # generate prime_list by reading all true
    prime_list = prime_bools.nonzero()[0]

# returns n divided by its smallest prime divisor
# TODO: is there a "cached" list pattern?
def get_reduction(n):

    global reductions

    if reductions[n]:
        return reductions[n]

    if prime_bools[n]:
        reductions[n] = n
        return n

    upper_bound = math.sqrt(abs(n))

    for i in np.nditer(prime_list):
        if i > upper_bound:
            return int(n)
        if n % i == 0:
            # TODO: single slash division when dividing?
            # if is divisible, do we still get float?
            return int(n // i)

    raise Exception("Prime table too small for input")

def gcd(a, b):
    while a:
        a, b = b % a, a
    return b

# handles the section before the first node
def generate_sequence(a, b):

    while True:

        if a == b and gcd(a, b) != 1:
            return None
        elif a % 2 and b % 2 and gcd(a, b) == 1:
            return generate_sequence_main(a, b)

        a, b = b, get_reduction(a + b)

# creates all requisite nodes and returns beginning
def generate_sequence_main(a, b):

    global cycles, nodes

    orig_a = a
    orig_b = b

    if (a,b) not in nodes:

        prev_prev = 0
        prev = 0

        # add to cache
        node = Node()
        nodes[(a,b)] = node

        # we enter node mode when preceded by odd and even
        while not (prev_prev % 2 and not prev % 2):

            node.terms += [a]
            node.length += 1
            prev_prev, prev, a, b = prev, a, b, get_reduction(a + b)

        # indicate succeeding node
        node.next = (a, b)

        # make sure succeeding nodes are created, and propagate end condition
        node_next = generate_sequence_main(a, b)
        node.cycle = node_next.cycle

    elif not nodes[(a,b)].cycle:

        # we have visited this node before
        # BUT we have not accounted for the corresponding cycle or we have not
        # therefore, construct the corresponding cycle and put it in the records
        cycle = Cycle()
        node = nodes[(a,b)]
        cycle.nodes.add((a,b))
        cycle.length += node.length

        while node.next != (a,b):
            cycle.nodes.add(node.next)
            node = nodes[node.next]
            cycle.length += node.length

        # naming convention
        name = str(cycle.length) + "-" + str(sorted(cycle.nodes)[0])

        cycles[name] = cycle
        nodes[(a,b)].cycle = name

    return nodes[(orig_a, orig_b)]


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('integers', metavar='N', type=int, nargs='+')
    args = parser.parse_args()

    generate_prime_lists(MAX_PRIME)

    if len(args.integers) == 1:

        # dictionary indexed by cycles
        cycle_counts = {"trivial": 0}
        bound = args.integers[0] + 1

        for a in range(1, bound):
            for b in range(1, bound):
                node = generate_sequence(a, b)
                if not node:
                    cycle_counts["trivial"] += 1
                else:
                    if node.cycle not in cycle_counts:
                        cycle_counts[node.cycle] = 0
                    cycle_counts[node.cycle] += 1

        print cycle_counts

    elif len(args.integers) == 2:

        a, b = tuple(args.integers)

        # if two terms: then just print one sequence
        node = generate_sequence(a, b)

        if not node:
            print "trivial"
            return

        # # print the node system
        # for index, cur_node in nodes.iteritems():
        #     print index, vars(cur_node)

        # # print the cycle list
        # for index, cycle in cycles.iteritems():
        #     print index, vars(cycle)

        print node.cycle

        # now construct output by traversal
        output = []
        visited = {}
        node_index = node.index()
        while node_index not in visited:
            node = nodes[node_index]
            output += node.terms
            visited[node_index] = True
            node_index = node.next

        # repeating is done to clarify the loop
        output += nodes[node_index].terms[:2] + ["..."]

        print ", ".join([str(x) for x in output])
        

if __name__ == "__main__":
    main()
