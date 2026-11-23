import argparse
import itertools
import json
import math
from pathlib import Path


ELEMENTS = tuple(itertools.product(range(4), range(2)))
IDENTITY = (0, 0)


def multiply(left, right):
    return ((left[0] + (-1) ** left[1] * right[0] + 2 * left[1] * right[1]) % 4,
            (left[1] + right[1]) % 2)


def inverse(element):
    return next(candidate for candidate in ELEMENTS
                if multiply(element, candidate) == IDENTITY
                and multiply(candidate, element) == IDENTITY)


def conjugate(conjugator, element):
    return multiply(multiply(conjugator, element), inverse(conjugator))


def generated_subgroup(subset):
    generated = {IDENTITY, *subset}
    while True:
        expanded = generated | {multiply(left, right) for left in generated for right in generated}
        if expanded == generated:
            return tuple(sorted(generated))
        generated = expanded


def is_good(subset, mapping):
    if any(mapping[mapping[element]] != element for element in subset):
        return False
    return all(mapping[conjugate(conjugator, element)] == conjugate(conjugator, mapping[element])
               and conjugate(mapping[conjugator], element) == conjugate(inverse(conjugator), element)
               for conjugator in subset for element in subset)


def analyse_subset(subset):
    subgroup = generated_subgroup(subset)
    center = tuple(element for element in subgroup
                   if all(multiply(element, other) == multiply(other, element) for other in subgroup))
    unseen = set(subset)
    classes = []
    while unseen:
        representative = min(unseen)
        orbit = frozenset(conjugate(element, representative) for element in subgroup)
        assert orbit <= unseen
        classes.append(orbit)
        unseen -= orbit
    class_index = {element: index for index, orbit in enumerate(classes) for element in orbit}
    destinations = []
    for orbit in classes:
        row = []
        for label in center:
            target_elements = {multiply(label, inverse(element)) for element in orbit}
            matches = [index for index, target in enumerate(classes) if target_elements == target]
            assert len(matches) <= 1
            row.append(matches[0] if matches else None)
        destinations.append(row)
    weights = [[row.count(index) for index in range(len(classes))] for row in destinations]
    assert all(weights[left][right] == weights[right][left]
               for left in range(len(classes)) for right in range(len(classes)))
    accepted = 0
    for labels in itertools.product(range(len(center)), repeat=len(classes)):
        targets = [destinations[index][label] for index, label in enumerate(labels)]
        if any(target is None or labels[target] != labels[index] for index, target in enumerate(targets)):
            continue
        mapping = {element: multiply(center[labels[class_index[element]]], inverse(element))
                   for element in subset}
        assert set(mapping.values()) <= set(subset) and is_good(subset, mapping)
        accepted += 1

    def weighted_count(remaining):
        if not remaining:
            return 1
        first, *rest = remaining
        return (weights[first][first] * weighted_count(tuple(rest))
                + sum(weights[first][partner]
                      * weighted_count(tuple(index for index in rest if index != partner))
                      for partner in rest))

    assert accepted == weighted_count(tuple(range(len(classes))))
    if len(center) >= 2:
        equality_criterion = all(target == index for index, row in enumerate(destinations) for target in row)
    else:
        equality_criterion = {inverse(element) for element in subset} == set(subset)
    assert (accepted == len(center) ** len(classes)) == equality_criterion
    direct = None
    if len(subset) <= 4:
        direct = sum(is_good(subset, dict(zip(subset, values)))
                     for values in itertools.product(subset, repeat=len(subset)))
        assert direct == accepted
    return {'elements': [list(element) for element in subset], 'generated_subgroup_order': len(subgroup),
            'center_order': len(center), 'orbit_count': len(classes), 'weights': weights,
            'good_involutions': accepted, 'weighted_count': accepted, 'direct_function_count': direct,
            'central_label_bound': len(center) ** len(classes), 'equality_criterion': equality_criterion}


def quaternion_example():
    subset = ((1, 0), (0, 1), (3, 0), (2, 1))
    for left, right in itertools.product(range(4), repeat=2):
        assert conjugate(subset[left], subset[right]) == subset[(2 * left - right) % 4]
    good = []
    for values in itertools.product(range(4), repeat=4):
        if is_good(subset, {element: subset[values[index]] for index, element in enumerate(subset)}):
            good.append(values)
    automorphisms = []
    antiautomorphisms = []
    for values in itertools.permutations(range(4)):
        mapping = {element: subset[values[index]] for index, element in enumerate(subset)}
        if all(mapping[conjugate(left, right)] == conjugate(mapping[left], mapping[right])
               for left in subset for right in subset):
            automorphisms.append(values)
        if all(mapping[conjugate(left, right)] == conjugate(inverse(mapping[left]), mapping[right])
               for left in subset for right in subset):
            antiautomorphisms.append(values)
    assert len(good) == 4 and len(automorphisms) == 8 and automorphisms == antiautomorphisms
    unseen = set(good)
    isomorphism_classes = []
    while unseen:
        representative = min(unseen)
        orbit = set()
        for automorphism in automorphisms:
            inverse_permutation = tuple(automorphism.index(index) for index in range(4))
            orbit.add(tuple(automorphism[representative[inverse_permutation[index]]] for index in range(4)))
        assert orbit <= set(good)
        unseen -= orbit
        isomorphism_classes.append([list(values) for values in sorted(orbit)])
    assert sorted(map(len, isomorphism_classes)) == [1, 1, 2]
    analysis = analyse_subset(tuple(sorted(subset)))
    assert (analysis['generated_subgroup_order'], analysis['center_order'], analysis['orbit_count']) == (8, 2, 2)
    assert min(len(automorphisms), analysis['central_label_bound']) == 4
    return {'coordinate_order': [list(element) for element in subset], 'dihedral_R4_identification': True,
            'all_functions_examined': 4 ** 4, 'good_maps': [list(values) for values in good],
            'automorphisms': [list(values) for values in automorphisms],
            'antiautomorphisms': [list(values) for values in antiautomorphisms],
            'symmetric_isomorphism_classes': isomorphism_classes, 'subquandle': analysis}


def linear_family():
    results = []
    for parameter in range(1, 32):
        modulus = 4 * parameter
        action = 2 * parameter - 1
        elements = tuple(range(modulus))
        assert math.gcd(action, modulus) == 1 and action ** 2 % modulus == 1
        operation = [[(action * target + (1 - action) * conjugator) % modulus
                      for target in elements] for conjugator in elements]
        fixed_points = tuple(element for element in elements if action * element % modulus == element)
        unseen = set(elements)
        orbits = []
        while unseen:
            orbit = {min(unseen)}
            frontier = list(orbit)
            while frontier:
                target = frontier.pop()
                for conjugator in elements:
                    image = operation[conjugator][target]
                    if image not in orbit:
                        orbit.add(image)
                        frontier.append(image)
            assert orbit <= unseen
            unseen -= orbit
            orbits.append(orbit)
        orbit_index = {element: index for index, orbit in enumerate(orbits) for element in orbit}
        accepted = []
        for labels in itertools.product(fixed_points, repeat=len(orbits)):
            mapping = tuple((labels[orbit_index[element]] - action * element) % modulus
                            for element in elements)
            if any(mapping[mapping[element]] != element for element in elements):
                continue
            assert all(mapping[operation[conjugator][target]] == operation[conjugator][mapping[target]]
                       and operation[mapping[conjugator]][target] == operation[conjugator][target]
                       for conjugator in elements for target in elements)
            accepted.append(mapping)
        assert len(accepted) == (10 if parameter % 2 else 4)
        if parameter % 2:
            assert len(orbits) == 4 and all(len(orbit) == parameter for orbit in orbits)
            component_permutations = set()
            for mapping in accepted:
                assert all(mapping[element] % parameter == element % parameter for element in elements)
                induced = tuple(mapping[next(element for element in elements if element % 4 == residue)] % 4
                                for residue in range(4))
                assert all(mapping[element] % 4 == induced[element % 4] for element in elements)
                assert all(induced[induced[residue]] == residue for residue in range(4))
                component_permutations.add(induced)
            assert len(component_permutations) == 10
        else:
            assert len(orbits) == 2 and len(fixed_points) == 2
        results.append({'parameter': parameter, 'order': modulus, 'action': action,
                        'fixed_points': list(fixed_points), 'component_count': len(orbits),
                        'component_sizes': [len(orbit) for orbit in orbits],
                        'label_assignments_examined': len(fixed_points) ** len(orbits),
                        'good_involutions': len(accepted)})
    return results


def main():
    parser = argparse.ArgumentParser(description='Exact quaternion and inverse-class checks.')
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    for element in ELEMENTS:
        assert multiply(element, IDENTITY) == element == multiply(IDENTITY, element)
        assert multiply(element, inverse(element)) == IDENTITY
    for left, middle, right in itertools.product(ELEMENTS, repeat=3):
        assert multiply(multiply(left, middle), right) == multiply(left, multiply(middle, right))
    analyses = []
    for size in range(1, len(ELEMENTS) + 1):
        for subset in itertools.combinations(ELEMENTS, size):
            if all(conjugate(left, right) in subset for left in subset for right in subset):
                analyses.append(analyse_subset(subset))
    whole = next(row for row in analyses if len(row['elements']) == 8)
    center_order = whole['center_order']
    involution_count = sum(math.factorial(center_order)
                           // (2 ** pairs * math.factorial(pairs) * math.factorial(center_order - 2 * pairs))
                           for pairs in range(center_order // 2 + 1))
    sharper_bound = involution_count * center_order ** (whole['orbit_count'] - center_order)
    assert whole['good_involutions'] == sharper_bound == 16 < whole['central_label_bound'] == 32
    first_class = sorted({pow(5, exponent, 16) for exponent in range(4)})
    inverse_class = sorted({(-exponent) % 16 for exponent in first_class})
    assert first_class == [1, 5, 9, 13] and inverse_class == [3, 7, 11, 15]
    result = {'schema_version': 1, 'quaternion_example': quaternion_example(),
              'whole_Q8': {'good_involutions': 16, 'sharper_bound': sharper_bound, 'central_label_bound': 32},
              'G5_inverse_classes': {'class_of_a': first_class, 'class_of_inverse_a': inverse_class},
              'nonempty_conjugation_closed_subsets_of_Q8': len(analyses), 'subset_checks': analyses,
              'linear_quandle_family': linear_family(),
              'coverage': 'Direct quaternion maps, automorphisms, antiautomorphisms and symmetric isomorphism classes; weighted counts and equality criterion on every quaternion conjugation subquandle; distinct inverse classes in G5; complete Alexander-label enumeration and defining identities for the linear family with parameters 1 through 31.'}
    output = json.dumps(result, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
