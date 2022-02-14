import loompy
from collections import namedtuple


class LoomReader:
    def __init__(self, filename):
        self.ds = loompy.connect(filename)

        self.cells = list(self.ds.ca['CellID'])
        self.genes = list(self.ds.ra['Gene'])

    def data(self):
        '''
        Returns all data for each cell.
        '''
        cells = {}
        data_zip = zip(
            self.ds.ca['Valid'],
            self.ds.ca['CellID'],
            self.ds.ca['_tSNE_1'],
            self.ds.ca['_tSNE_2'],
            self.ds.ca['ClusterID'],
            self.ds.ca['ClusterName'],
            self.ds.ca['X'],
            self.ds.ca['Y']
        )
        for (valid, cell_id, tsne1, tsne2,
             cluster_id, cluster_name, x, y) in data_zip:
            if valid:
                cells[cell_id] = {
                    'mappings': {
                        't-SNE': [tsne1, tsne2]
                    },
                    'cluster': cluster_name,
                    'genes': self.by_cell(cell_id),
                    'xy': [x, y]
                }
        return cells

    def clusters(self):
        '''
        Returns information about each cluster.

        >>> lr = LoomReader('fake-files/input/linnarsson/\
linnarsson.cells.loom')
        >>> clusters = lr.clusters()
        >>> clusters[4].name
        'pyramidal L4'
        >>> sorted(clusters[4].cell_ids)[:3]
        ['104', '110', '111']
        '''
        clusters = {}
        Cluster = namedtuple('Cluster', ['name', 'cell_ids'])

        clustered_cell_zip = zip(
            self.ds.ca['Valid'],
            self.ds.ca['ClusterID'],
            self.ds.ca['ClusterName'],
            self.ds.ca['CellID']
        )
        for (valid, cluster_id, cluster_name, cell_id) in clustered_cell_zip:
            if valid:
                if cluster_id in clusters:
                    clusters[cluster_id].cell_ids.append(cell_id)
                else:
                    clusters[cluster_id] = Cluster(cluster_name, [cell_id])
        return clusters

    def tsne(self):
        '''
        Returns a tSNE pair for each valid cell.

        >>> lr = LoomReader('fake-files/input/linnarsson/\
linnarsson.cells.loom')
        >>> tsne = lr.tsne()
        >>> tsne['42']
        (-6.917..., 16.644...)
        '''
        cells = {}

        tsne_zip = zip(
            self.ds.ca['Valid'],
            self.ds.ca['CellID'],
            self.ds.ca['_tSNE_1'],
            self.ds.ca['_tSNE_2']
        )
        for (valid, cell_id, tsne1, tsne2) in tsne_zip:
            if valid:
                cells[cell_id] = (tsne1, tsne2)
        return cells

    def xy(self):
        '''
        Returns the xy position for each valid cell.

        >>> lr = LoomReader('fake-files/input/linnarsson/\
linnarsson.cells.loom')
        >>> xy = lr.xy()
        >>> xy['42']
        (21164.7..., 35788.1...)
        '''
        cells = {}

        tsne_zip = zip(
            self.ds.ca['Valid'],
            self.ds.ca['CellID'],
            self.ds.ca['X'],
            self.ds.ca['Y']
        )
        for (valid, cell_id, x, y) in tsne_zip:
            if valid:
                cells[cell_id] = (x, y)
        return cells

    def by_cell(self, cell_id):
        '''
        Given a cell_id, returns a dict with the values for each gene.

        >>> lr = LoomReader('fake-files/input/linnarsson/\
linnarsson.cells.loom')
        >>> lr.by_cell('42')['Gad2']
        7
        '''
        return dict(zip(
            self.genes,
            (int(x[0]) for x in self.ds[:, self.ds.ca.CellID == cell_id])
        ))

    def by_gene(self, gene):
        '''
        Given a gene, returns a dict with the values for each cell.

        >>> lr = LoomReader('fake-files/input/linnarsson/\
linnarsson.cells.loom')
        >>> lr.by_gene('Gad2')['42']
        7
        '''
        return dict(zip(
            self.cells,
            self.ds[self.ds.ra.Gene == gene, :].tolist()[0]
        ))
