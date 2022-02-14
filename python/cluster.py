from scipy.cluster.hierarchy import leaves_list, linkage
import pandas as pd

# I found this tutorial very helpful:
# https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/


def _order_rows(dataframe):
    '''
    >>> df = pd.DataFrame({
    ...   'cell-1': {'a':8, 'b':1, 'c': 7, 'd': 2},
    ...   'cell-2': {'a':1, 'b':1, 'c': 1, 'd': 1},
    ...   'cell-3': {'a':9, 'b':1, 'c': 8, 'd': 2},
    ...   'cell-4': {'a':1, 'b':2, 'c': 1, 'd': 1}
    ... })
    >>> _order_rows(df)
    ['a', 'c', 'b', 'd']

    '''
    row_labels = dataframe.index.tolist()
    if len(row_labels) > 1:
        rows_linkage = linkage(dataframe, 'ward')
        rows_order = leaves_list(rows_linkage).tolist()
        return [row_labels[i] for i in rows_order]
    else:
        return row_labels


def _order(dataframe):
    '''
    >>> df = pd.DataFrame({
    ...   'cell-1': {'a':8, 'b':1, 'c': 7},
    ...   'cell-2': {'a':1, 'b':1, 'c': 1},
    ...   'cell-3': {'a':9, 'b':1, 'c': 8}
    ... })
    >>> _order(df)['rows']
    ['b', 'a', 'c']
    >>> _order(df)['cols']
    ['cell-2', 'cell-1', 'cell-3']

    '''
    col_label_order = _order_rows(dataframe.T)
    row_label_order = _order_rows(dataframe)
    return {'rows': row_label_order, 'cols': col_label_order}


def _to_dataframe(cells):
    '''
    >>> cells = {
    ...   'cell-1': { 'genes': {'a':8, 'b':1, 'c': 7}, 'extra': 'f'},
    ...   'cell-2': { 'genes': {'a':1, 'b':1, 'c': 1}, 'extra': 'field'},
    ...   'cell-3': { 'genes': {'a':9, 'b':1, 'c': 10}, 'extra': 'field'}
    ... }
    >>> _to_dataframe(cells)
       cell-1  cell-2  cell-3
    a     0.8     0.1     0.9
    b     0.1     0.1     0.1
    c     0.7     0.1     1.0

    '''
    clean = {}
    for k, v in cells.items():
        clean[k] = v['genes']
    df = pd.DataFrame(clean)
    df_max = df.values.max()
    # If we don't round, small differences in the floating point representation
    # cause cluster results on Travis to be slightly different than local.
    return (df / df_max).round(4)


def _row_norm(df):
    '''
    >>> df = pd.DataFrame({
    ...   'cell-1': {'a':8, 'b':1, 'c': 7},
    ...   'cell-2': {'a':1, 'b':1, 'c': 1},
    ...   'cell-3': {'a':10, 'b':1, 'c': 10}
    ... })
    >>> _row_norm(df)
       cell-1  cell-2  cell-3
    a     0.8     0.1     1.0
    b     1.0     1.0     1.0
    c     0.7     0.1     1.0
    '''
    t = df.T
    return (t / t.max()).T


def cluster(cells):
    '''
    >>> cells = {
    ...   'cell-1': { 'genes': {'a':8, 'b':2, 'c': 7}, 'extra': 'field'},
    ...   'cell-2': { 'genes': {'a':1, 'b':1, 'c': 1}, 'extra': 'field'},
    ...   'cell-3': { 'genes': {'a':10, 'b':2, 'c': 10}, 'extra': 'field'}
    ... }
    >>> clustered = cluster(cells)
    >>> clustered['rows']
    ['b', 'a', 'c']
    >>> clustered['cols']
    ['cell-2', 'cell-1', 'cell-3']
    >>> clustered['matrix']
    [[0.5, 1.0, 1.0], [0.1, 0.8, 1.0], [0.1, 0.7, 1.0]]

    '''
    df = _to_dataframe(cells)
    rows_cols = _order(df)
    clustered = df[rows_cols['cols']].loc[rows_cols['rows']]
    return {
        'rows': rows_cols['rows'],
        'cols': rows_cols['cols'],
        'matrix': _row_norm(clustered).round(3).values.tolist()
    }
