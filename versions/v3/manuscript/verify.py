import argparse
import itertools
import json
import math
from pathlib import Path


def analyse(action):
    elements = tuple(itertools.product(range(16), range(4)))
    identity = (0, 0)

    def multiply(left, right):
        return ((left[0] + pow(action, left[1], 16) * right[0]) % 16,
                (left[1] + right[1]) % 4)

    def inverse(element):
        return ((-pow(action, -element[1], 16) * element[0]) % 16,
                -element[1] % 4)

    def conjugate(conjugator, element):
        return multiply(multiply(conjugator, element), inverse(conjugator))

    assert pow(action, 4, 16) == 1
    for element in elements:
        assert multiply(element, inverse(element)) == identity
        assert multiply(inverse(element), element) == identity
        assert multiply(element, identity) == element == multiply(identity, element)
    for left, middle, right in itertools.product(elements, repeat=3):
        assert multiply(multiply(left, middle), right) == multiply(left, multiply(middle, right))
    assert multiply((1, 0), (0, 1)) != multiply((0, 1), (1, 0))
    center = tuple(element for element in elements
                   if all(multiply(element, other) == multiply(other, element) for other in elements))
    classes = []
    unseen = set(elements)
    while unseen:
        representative = min(unseen)
        conjugacy_class = frozenset(conjugate(conjugator, representative) for conjugator in elements)
        assert conjugacy_class <= unseen
        classes.append(conjugacy_class)
        unseen -= conjugacy_class
    class_index = {element: index for index, conjugacy_class in enumerate(classes)
                   for element in conjugacy_class}
    for conjugator, element in itertools.product(elements, repeat=2):
        coordinate_formula = (((1 - pow(action, element[1], 16)) * conjugator[0]
                               + pow(action, conjugator[1], 16) * element[0]) % 16,
                              element[1])
        assert conjugate(conjugator, element) == coordinate_formula
    destinations = []
    for conjugacy_class in classes:
        row = []
        for central in center:
            targets = {class_index[multiply(central, inverse(element))] for element in conjugacy_class}
            assert len(targets) == 1
            row.append(targets.pop())
        destinations.append(row)
    for index, row in enumerate(destinations):
        for label, target in enumerate(row):
            assert destinations[target][label] == index
    unseen_classes = set(range(len(classes)))
    blocks = []
    while unseen_classes:
        component = {min(unseen_classes)}
        frontier = list(component)
        while frontier:
            for target in destinations[frontier.pop()]:
                if target not in component:
                    component.add(target)
                    frontier.append(target)
        indices = sorted(component)
        unseen_classes -= component
        local_index = {index: position for position, index in enumerate(indices)}
        matrix = [[destinations[index].count(target) for target in indices] for index in indices]
        assert all(matrix[row][column] == matrix[column][row]
                   for row, column in itertools.product(range(len(indices)), repeat=2))
        accepted = 0
        for labels in itertools.product(range(len(center)), repeat=len(indices)):
            if any(labels[local_index[destinations[index][labels[position]]]] != labels[position]
                   for position, index in enumerate(indices)):
                continue
            accepted += 1
            domain = set().union(*(classes[index] for index in indices))
            rho = {element: multiply(center[labels[local_index[class_index[element]]]], inverse(element))
                   for element in domain}
            assert set(rho.values()) == domain
            assert all(rho[rho[element]] == element for element in domain)
            for conjugator, element in itertools.product(elements, domain):
                assert rho[conjugate(conjugator, element)] == conjugate(conjugator, rho[element])
            for element, test in itertools.product(domain, elements):
                assert conjugate(rho[element], test) == conjugate(inverse(element), test)

        def matching_count(remaining):
            if not remaining:
                return 1
            first, *rest = remaining
            total = matrix[first][first] * matching_count(tuple(rest))
            for partner in rest:
                total += matrix[first][partner] * matching_count(
                    tuple(index for index in rest if index != partner))
            return total

        assert accepted == matching_count(tuple(range(len(indices))))
        blocks.append({'representatives': [list(min(classes[index])) for index in indices],
                       'class_sizes': [len(classes[index]) for index in indices],
                       'weights': matrix, 'count': accepted})
    return {'action': action, 'order': len(elements), 'center': [list(element) for element in center],
            'center_order': len(center), 'class_number': len(classes),
            'blocks': blocks, 'good_involutions': math.prod(block['count'] for block in blocks)}


def main():
    parser = argparse.ArgumentParser(description='Exact good-involution counts for two semidirect products.')
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    results = [analyse(action) for action in (5, 7)]
    assert [(entry['order'], entry['center_order'], entry['class_number']) for entry in results] == [(64, 4, 22)] * 2
    assert [entry['good_involutions'] for entry in results] == [8847360, 2160000]
    output = json.dumps({'schema_version': 1, 'groups': results,
                         'coverage': 'Exact group laws, centers, conjugacy classes, central-label blocks, and complete good-involution counts.'}, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
