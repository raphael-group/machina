#!/usr/env python

def adjusted_rand_index_pair_counts(a, b, c, d):
    """
    Compute the adjusted Rand index from pair counts; helper function.

    Arguments:
    a: number of pairs of elements that are clustered in both partitions
    b: number of pairs of elements that are clustered in first but not second partition
    c: number of pairs of elements that are clustered in second but not first partition
    d: number of pairs of elements that are clustered in neither partiition

    Example usage:
    In [1]: a = 1
    In [2]: b = 2
    In [3]: c = 2
    In [4]: d = 10
    In [5]: adjusted rand_index_pair_counts(a, b, c, d)
    Out[5]: 0.16666666666666666
    """
    if a*(b+c+2*d)+b*(b+d)+c*(c+d)!=0:
        return float(2*(a*d-b*c))/float(a*(b+c+2*d)+b*(b+d)+c*(c+d))
    else:
        return 1.0

def recall_and_precision(X, Y):
    """
    Compute the recall and precision from partitions.

    Arguments:
    X: inferred partition given as list of lists
    Y: true partition given as list of lists

    Example usage:
    In [1]: U = [[0, 1], [2, 3], [4, 5]]
    In [2]: V = [[0], [1], [2], [3, 4, 5]]
    In [3]: recall_and_precision(U, V)
    Out[3]: (0.09090909090909091, 0.3333333333333333)
    """
    n = sum(len(x) for x in X)

    TP, FP, FN, TN = count_pairs(X, Y, n)

    recall = float(TP) / float(TP + FN)
    precision = float(TP) / float(TP + FP)

    return recall, precision

def adjusted_rand_index(X, Y):
    """
    Compute the adjusted Rand index from partitions.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists

    Example usage:
    In [1]: U = [[0, 1], [2, 3], [4, 5]]
    In [2]: V = [[0], [1], [2], [3, 4, 5]]
    In [3]: adjusted_rand_index(U, V)
    Out[3]: 0.16666666666666666
    """
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)
    return adjusted_rand_index_pair_counts(a, b, c, d)

def cluster_similarity_scores(X, Y, measures):
    """
    Compute several cluster similarity scores.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists
    measures: measures given as list of strings

    Example usage:
    In [1]: X = [[0, 1], [2, 3], [4, 5]]
    In [2]: Y = [[0], [1], [2], [3, 4, 5]]
    In [3]: count_pairs(X, Y, ['ri', 'ari'])
    Out[3]: [0.7333333333333333, 0.16666666666666666]
    """

    # Check measures.
    implemented_measures = set(['ji', 'ri', 'ari', 'fmi', 'f1'])
    for measure in measures:
        if measure not in implemented_measures:
            raise NotImplementedError('{} has not been implemented; choose from {}.'.format(measure, ', '.join(sorted(implemented_measures))))

    # Compute pair counts.
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)

    # Compute cluster similarity scores.
    results = [0.0 for measure in measures]
    for k, measure in enumerate(measures):
        if measure=='ji':
            results[k] = jaccard_index_pair_counts(a, b, c)
        elif measure=='ri':
            results[k] = rand_index_pair_counts(a, b, c, d)
        elif measure=='ari':
            results[k] = adjusted_rand_index_pair_counts(a, b, c, d)
        elif measure=='fmi':
            results[k] = fowlkes_mallows_index_pair_counts(a, b, c)
        elif measure=='f1':
            results[k] = f1_score_pair_counts(a, b, c)
    return results

def count_pairs(X, Y, n):
    """
    Find the number of pairs of elements that belong to similar or different clusters in different
    partitions.

    Arguments:
    U: partition of elements given as a list of lists
    V: partition of elements given as a list of lists
    n: number of elements

    Example usage:
    In [1]: X = [[0, 1], [2, 3], [4, 5]]
    In [2]: Y = [[0], [1], [2], [3, 4, 5]]
    In [3]: n = 6
    In [4]: count_pairs(X, Y, n)
    Out[4]: (1, 10, 2, 2)
    """
    import numpy as np

    # We find pairs of element that belong to similar or different clusters in different
    # partitions. For a pair of elements to belong to the same cluster, the cluster must have at
    # least two elements, so it is enough to consider pairs of elements belonging to nonsingleton
    # clusters.
    U = [set(x) for x in X if len(x)>1]
    V = [set(y) for y in Y if len(y)>1]

    # Find the number of pairs of elements that are clustered together in the first partition.
    P = np.array([len(u) for u in U])
    p = int(np.sum(P*(P-1))/2)

    # Find the number of pairs of elements that are clustered together in the second partition.
    Q = np.array([len(v) for v in V])
    q = int(np.sum(Q*(Q-1))/2)

    # Find the number of pairs of elements that are clustered together in both partitions.
    R = np.array([len(set.intersection(u, v)) for v in V for u in U])
    r = int(np.sum(R*(R-1))/2)

    # Find the number of pairs of elements that are clustered together (a) in both partitions, (b)
    # the first partition but not the second partition, (c) the second partition but not the first
    # partition, and (d) neither partition.
    a = r
    b = p-a
    c = q-a
    d = int(n*(n-1)/2)-a-b-c

    return a, b, c, d

def f1_score_pair_counts(a, b, c):
    """
    Compute the F1 score from pair counts; helper function.

    Arguments:
    a: number of pairs of elements that are clustered in both partitions
    b: number of pairs of elements that are clustered in first but not second partition
    c: number of pairs of elements that are clustered in second but not first partition

    Example usage:
    In [1]: a = 1
    In [2]: b = 2
    In [3]: c = 2
    In [4]: d = 10
    In [5]: f1_score_pair_count(a, b, c)
    Out[5]: 0.3333333333333333
    """
    if 2*a+b+c!=0:
        return float(2*a)/float(2*a+b+c)
    else:
        return 1.0

def f1_score(X, Y):
    """
    Compute the F1 score from partitions.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists

    Example usage:
    In [1]: U = [[0, 1], [2, 3], [4, 5]]
    In [2]: V = [[0], [1], [2], [3, 4, 5]]
    In [3]: f1_score(U, V)
    Out[3]: 0.3333333333333333
    """
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)
    return f1_score_pair_counts(a, b, c)

def fowlkes_mallows_index_pair_counts(a, b, c):
    """
    Compute the Fowlkes-Mallows index from pair counts; helper function.

    Arguments:
    a: number of pairs of elements that are clustered in both partitions
    b: number of pairs of elements that are clustered in first but not second partition
    c: number of pairs of elements that are clustered in second but not first partition

    Example usage:
    In [1]: a = 1
    In [2]: b = 2
    In [3]: c = 2
    In [4]: d = 10
    In [5]: fowlkes_mallows_index_pair_count(a, b, c)
    Out[5]: 0.3333333333333333
    """
    import math

    if (a+b)*(a+c)!=0:
        return float(a)/math.sqrt((a+b)*(a+c))
    elif a+b==a+c==0:
        return 1.0
    else:
        return 0.0

def fowlkes_mallows_index(X, Y):
    """
    Compute the Fowlkes-Mallows index from partitions.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists

    Example usage:
    In [1]: X = [[0, 1], [2, 3], [4, 5]]
    In [2]: Y = [[0], [1], [2], [3, 4, 5]]
    In [3]: fowlkes_mallows_index(X, Y)
    Out[3]: 0.3333333333333333
    """
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)
    return fowlkes_mallows_index_pair_counts(a, b, c)

def jaccard_index_pair_counts(a, b, c):
    """
    Compute the Jaccard index from pair counts; helper function.

    Arguments:
    a: number of pairs of elements that are clustered in both partitions
    b: number of pairs of elements that are clustered in first but not second partition
    c: number of pairs of elements that are clustered in second but not first partition
    d: number of pairs of elements that are clustered in neither partition

    Example usage:
    In [1]: a = 1
    In [2]: b = 2
    In [3]: c = 2
    In [4]: d = 10
    In [5]: jaccard_index_pair_counts(a, b, c, d)
    Out[5]: 0.2
    """
    if a+b+c!=0:
        return float(a)/float(a+b+c)
    else:
        return 1.0

def jaccard_index(X, Y):
    """
    Compute the Jaccard index from partitions.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists

    Example usage:
    In [1]: X = [[0, 1], [2, 3], [4, 5]]
    In [2]: Y = [[0], [1], [2], [3, 4, 5]]
    In [3]: jaccard_index(X, Y, n)
    Out[3]: 0.2
    """
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)
    return jaccard_index_pair_counts(a, b, c)

def labels_to_lists(X):
    """
    Convert element labels to list of lists of elements.

    Argument:
    X: element labels

    Example usage:
    In [1]: labels_to_lists([3, 1, 4, 1, 5])
    Out[1]: [[1, 3], [0], [2], [4]]
    """
    clusters = set(X)
    index = {x: i for i, x in enumerate(clusters)}
    Y = [[] for _ in xrange(len(clusters))]
    for i, x in enumerate(X):
        Y[index[x]].append(i)
    return Y

def lists_to_labels(X):
    """
    Convert list of lists of elements to element labels.

    Argument:
    X: list of lists of elements

    Example usage:
    In [1]: lists_to_labels([[1, 5], [3], [2, 0, 4]])
    Out[1]: [2, 0, 2, 1, 2, 0]
    """
    elements = set(a for b in X for a in b)
    index = {a: i for i, a in enumerate(elements)}
    Y = [-1 for _ in xrange(len(elements))]
    for i, x in enumerate(X):
        for j in x:
            Y[index[j]] = i
    return Y

def rand_index_pair_counts(a, b, c, d):
    """
    Compute the Rand index from pair counts; helper function.

    Arguments:
    a: number of pairs of elements that are clustered in both partitions
    b: number of pairs of elements that are clustered in first but not second partition
    c: number of pairs of elements that are clustered in second but not first partition
    d: number of pairs of elements that are clustered in neither partiition

    Example usage:
    In [1]: a = 1
    In [2]: b = 2
    In [3]: c = 2
    In [4]: d = 10
    In [5]: rand_index_pair_counts(a, b, c, d)
    Out[5]: 0.7333333333333333
    """
    if a+b+c+d!=0:
        return float(a+d)/float(a+b+c+d)
    else:
        return 1.0

def rand_index(X, Y):
    """
    Compute the Rand index from partitions.

    Arguments:
    X: partition given as list of lists
    Y: partition given as list of lists

    Example usage:
    In [1]: X = [[0, 1], [2, 3], [4, 5]]
    In [2]: Y = [[0], [1], [2], [3, 4, 5]]
    In [3]: rand_index(X, Y)
    Out[3]: 0.7333333333333333
    """
    n = sum(len(x) for x in X)
    a, b, c, d = count_pairs(X, Y, n)
    return rand_index_pair_counts(a, b, c, d)
