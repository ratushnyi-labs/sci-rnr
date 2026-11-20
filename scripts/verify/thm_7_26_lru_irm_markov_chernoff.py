#!/usr/bin/env python3
"""Check Theorem 7.26 using finite-state LRU models and exact path sums."""

from collections import defaultdict
from itertools import permutations, product
from math import exp, prod

import numpy as np


def advance(stack, request, capacity):
    return (request, *(item for item in stack if item != request))[:capacity]


def stationary_model(popularity, capacity):
    alphabet = range(len(popularity))
    stacks = tuple(permutations(alphabet, capacity))
    states = tuple(product(stacks, alphabet))
    index = {state: i for i, state in enumerate(states)}
    transition = np.zeros((len(states), len(states)))
    for i, (stack, request) in enumerate(states):
        next_stack = advance(stack, request, capacity)
        for next_request, probability in enumerate(popularity):
            transition[i, index[next_stack, next_request]] = probability
    stationarity = np.vstack((transition.T - np.eye(len(states)), np.ones(len(states))))
    stationary = np.linalg.lstsq(stationarity, np.r_[np.zeros(len(states)), 1], rcond=None)[0]
    np.testing.assert_allclose(stationary @ transition, stationary, atol=1e-12)
    assert min(stationary) > 0

    weights = np.sqrt(stationary)
    scaled = weights[:, None] * transition / weights[None, :]
    eigenvalues = np.linalg.eigvalsh((scaled + scaled.T) / 2)
    gap = 1.0 if len(states) == 1 else 1 - max(0.0, eigenvalues[-2])
    steps = capacity + 1
    minorization = min(popularity) ** steps
    assert np.linalg.matrix_power(transition, steps).min() >= minorization - 1e-12
    assert gap >= minorization / steps**2 - 1e-12

    stack_law = {stack: sum(stationary[index[stack, x]] for x in alphabet) for stack in stacks}
    for stack, request in states:
        assert abs(stationary[index[stack, request]] - stack_law[stack] * popularity[request]) < 1e-12
    miss_rate = sum(stationary[i] for i, (stack, request) in enumerate(states) if request not in stack)
    return gap, stack_law, miss_rate


def miss_distribution(popularity, capacity, stack_law, queries):
    law = {(stack, 0): probability for stack, probability in stack_law.items()}
    for _ in range(queries):
        next_law = defaultdict(float)
        for (stack, misses), weight in law.items():
            for request, probability in enumerate(popularity):
                next_stack = advance(stack, request, capacity)
                next_law[next_stack, misses + int(request not in stack)] += weight * probability
        law = next_law
    counts = [0.0] * (queries + 1)
    for (_, misses), probability in law.items():
        counts[misses] += probability
    assert abs(sum(counts) - 1.0) < 1e-12
    return counts


def outside(counts, queries, mean, threshold):
    return sum(weight for misses, weight in enumerate(counts)
               if abs(misses / queries - mean) + 1e-14 >= threshold)


def check(popularity, capacity):
    gap, stationary, miss_rate = stationary_model(popularity, capacity)
    initial = next(iter(stationary))
    partial = initial[:-1]
    density = 1 / stationary[initial]
    block_probability = prod(popularity[item] for item in initial)
    active = 0
    for queries in (1, 8, 32, 128):
        equilibrium = miss_distribution(popularity, capacity, stationary, queries)
        warm = miss_distribution(popularity, capacity, {initial: 1.0}, queries)
        cold = miss_distribution(popularity, capacity, {partial: 1.0}, queries)
        assert abs(sum(m * weight for m, weight in enumerate(equilibrium)) / queries - miss_rate) < 1e-12
        for step in range(1, 20):
            epsilon = step / 20
            bound = 2 * exp(-queries * gap * epsilon**2)
            assert outside(equilibrium, queries, miss_rate, epsilon) <= bound + 1e-12
            assert outside(warm, queries, miss_rate, epsilon) <= density * bound + 1e-12
            for burn_in in (0, min(2 * capacity, queries), queries // 4):
                cold_bound = bound + (1 - block_probability) ** (burn_in // capacity)
                tail = outside(cold, queries, miss_rate, epsilon + burn_in / queries)
                assert tail <= cold_bound + 1e-12
                active += int(0 < tail and cold_bound < 1)
    print(f"N={len(popularity)}, C={capacity}: gap={gap:.6f}, stationary miss rate={miss_rate:.6f}, "
          f"nontrivial cold-cache checks={active}")
    return active


if __name__ == "__main__":
    models = [((1.0,), 1), ((0.5, 0.5), 1), ((0.55, 0.28, 0.17), 1),
              ((0.55, 0.28, 0.17), 2), ((0.5, 0.25, 0.15, 0.10), 3)]
    assert sum(check(popularity, capacity) for popularity, capacity in models) > 0
    print("PASS: stationary, fixed full-stack, and initially partial/empty LRU tail bounds.")
