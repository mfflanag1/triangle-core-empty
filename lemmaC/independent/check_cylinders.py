"""Exact finite-cylinder check of the rerooting calculation in Theorem 6.7."""

from fractions import Fraction
from itertools import product
from math import prod


# These three whole teams form a connected piece of the incidence tree.
# Every agent has two named team choices leading out of the piece.
TEAMS = ((0, 1, 2), (0, 3, 4), (1, 5, 6))
AGENTS = tuple(range(7))


def incident_weights(team_weights):
    weights = []
    for v in AGENTS:
        row = {e: team_weights[e] for e, team in enumerate(TEAMS) if v in team}
        remaining = 1 - sum(row.values())
        share = Fraction(v + 1, v + 3)
        row[("exit", v, 0)] = remaining * share
        row[("exit", v, 1)] = remaining * (1 - share)
        weights.append(row)
    return weights


def orientations(root):
    """List the parent-to-child transitions when the chosen agent is the root."""
    transitions = []
    seen = {root}
    pending = [root]
    while pending:
        parent = pending.pop()
        for team, members in enumerate(TEAMS):
            if parent in members:
                for child in members:
                    if child not in seen:
                        transitions.append((parent, team, child))
                        seen.add(child)
                        pending.append(child)
    return transitions


def chain_mass(states, weights, root, transitions):
    mass = weights[root][states[root]]
    for parent, edge, child in transitions:
        if not mass:
            return mass
        if states[parent] == edge:
            if states[child] != edge:
                return Fraction(0)
        elif states[child] == edge:
            return Fraction(0)
        else:
            mass *= weights[child][states[child]] / (1 - weights[parent][edge])
    return mass


def cylinder_mass(states, weights, team_weights):
    mass = prod(weights[v][states[v]] for v in AGENTS)
    if not mass:
        return mass
    for e, members in enumerate(TEAMS):
        chosen = sum(states[v] == e for v in members)
        if chosen not in (0, 3):
            return Fraction(0)
        prob = team_weights[e] if chosen else 1 - team_weights[e]
        mass /= prob**2
    return mass


def verify(team_weights):
    weights = incident_weights(team_weights)
    rooted = [orientations(root) for root in AGENTS]
    marginals = [{state: Fraction(0) for state in row} for row in weights]
    total = Fraction(0)
    count = 0
    for states in product(*(row.keys() for row in weights)):
        expected = cylinder_mass(states, weights, team_weights)
        for root in AGENTS:
            actual = chain_mass(states, weights, root, rooted[root])
            if actual != expected:
                raise ArithmeticError((root, states, actual, expected))
        total += expected
        for v in AGENTS:
            marginals[v][states[v]] += expected
        count += 1
    if total != 1 or marginals != weights:
        raise ArithmeticError("failed normalization or one-agent marginals")
    return count


if __name__ == "__main__":
    interior = (Fraction(2, 7), Fraction(3, 7), Fraction(1, 5))
    boundary = (Fraction(0), Fraction(1), Fraction(1, 3))
    counts = [verify(params) for params in (interior, boundary)]
    print(f"Verified {sum(counts)} state patterns at all {len(AGENTS)} roots.")
    print("Both laws have unit mass and the specified exact marginals.")
