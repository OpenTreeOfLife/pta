import ivy, json
from ivy import treegraph as tg
gt = tg.gt

g = tg.load_taxonomy_graph('applications/pta/static/ott2.8.gt.gz')

def buildtree(t, otu_id2data):
    from ivy.tree import Node
    node_id2data = t['nodeById']
    root = None
    for i, d in node_id2data.iteritems():
        n = Node()
        n.snode_id = i
        n.taxid = None
        if d.get('@root'):
            n.isroot = True
            root = n
        oid = d.get('@otu')
        if oid:
            n.isleaf = True
            try:
                n.otu = otu_id2data[oid]
                n.label = (n.otu.get('^ot:ottTaxonName') or
                           n.otu.get('^ot:originalLabel'))
                n.taxid = n.otu.get('^ot:ottId')
            except KeyError:
                pass
                ## print t['nexson_file'], 'missing', oid
        n.nexson_id = i
        d['node'] = n
    for nid, ed in t['edgeBySourceId'].iteritems():
        n = node_id2data[nid]
        for e in ed.itervalues():
            cd = node_id2data[e['@target']]
            c = cd['node']
            c.length = e.get('@length')
            n['node'].add_child(c)
    root.tree_nexson = t
    root.stree = t['treeid']
    root.ladderize()
    ivy.tree.index(root)
    return root

def proctree(r, i=1):
    tg.map_stree(g, r)
    taxids = set()
    for lf in r.leaves():
        taxids.update(lf.taxid_rootpath)
    taxg = tg.taxid_new_subgraph(g, taxids)
    verts = taxg.new_vertex_property('bool')
    edges = taxg.new_edge_property('bool')

    # add stree's nodes and branches into taxonomy graph
    tg.merge_stree(taxg, r, i, verts, edges)

    # next, add taxonomy edges to taxg connecting 'incertae sedis'
    # leaves in stree to their containing taxa
    for lf in r.leaves():
        if lf.taxid and lf.taxid in taxg.taxid_vertex and lf.incertae_sedis:
            taxv = taxg.taxid_vertex[lf.taxid]
            ev = taxg.edge(taxv, lf.v, True)
            if ev:
                assert len(ev)==1
                e = ev[0]
            else:
                e = taxg.add_edge(taxv, lf.v)
            taxg.edge_in_taxonomy[e] = 1
    
    # make a view of taxg that keeps only the vertices and edges traced by
    # the source tree
    gv = tg.graph_view(taxg, vfilt=verts, efilt=edges)
    gv.vertex_strees = taxg.vertex_strees
    gv.edge_strees = taxg.edge_strees

    # the following code sets up the visualization
    ecolor = taxg.new_edge_property('string')
    for e in taxg.edges():
        est = taxg.edge_strees[e]
        eit = taxg.edge_in_taxonomy[e]
        if len(est) and not eit: ecolor[e] = 'blue'
        elif len(est) and eit: ecolor[e] = 'green'
        else: ecolor[e] = 'yellow'

    ewidth = taxg.new_edge_property('int')
    for e in taxg.edges():
        est = taxg.edge_strees[e]
        if len(est): ewidth[e] = 3
        else: ewidth[e] = 1

    vcolor = taxg.new_vertex_property('string')
    for v in taxg.vertices():
        if not taxg.vertex_in_taxonomy[v]: vcolor[v] = 'blue'
        else: vcolor[v] = 'green'

    vsize = taxg.new_vertex_property('int')
    for v in taxg.vertices():
        if taxg.vertex_in_taxonomy[v] or v.out_degree()==0:
            vsize[v] = 4
        else: vsize[v] = 2

    pos, pin = tg.layout(taxg, gv, gv.root, sfdp=True, deg0=195.0,
                         degspan=150.0, radius=400)

    for v in gv.vertices(): pin[v] = 1

    for e in taxg.edges():
        src = e.source()
        tgt = e.target()
        if not verts[src]:
            verts[src] = 1
            pos[src] = [0.0, 0.0]
            vcolor[src] = 'red'
        if not verts[tgt]:
            verts[tgt] = 1
            pos[tgt] = [0.0, 0.0]
            vcolor[tgt] = 'red'
        if not edges[e]:
            edges[e] = 1
            ecolor[e] = 'red'
            ewidth[e] = 1.0
            gv.wt[e] = 1.0

    pos = gt.sfdp_layout(gv, pos=pos, pin=pin, eweight=gv.wt, multilevel=False)

    nodes = []
    links = []
    idx = {}

    xmin = min([ pos[x][0] for x in gv.vertices() ])
    ymin = min([ pos[x][1] for x in gv.vertices() ])
    for x in gv.vertices(): pos[x] = [pos[x][0]-xmin, pos[x][1]-ymin]

    for i,v in enumerate(gv.vertices()):
        idx[int(v)] = i
        taxid = gv.vertex_taxid[v]
        name = gv.taxid_name(taxid) if taxid else 'node%s' % int(v)
        isleaf = v.out_degree()==0
        d = dict(label=name, isleaf=isleaf, strees=list(gv.vertex_strees[v]))
        if taxid: d['taxid'] = taxid
        ## if dist: d['dist'] = dist[v]
        if pos and pos[v]:
            x, y = pos[v]
            d['x'] = x; d['y'] = y
            d['fixed'] = True
        d['color'] = vcolor[v]
        d['size'] = vsize[v]
        ## if vtext: d['label'] = vtext[v]
        d['altlabel'] = gv.vertex_name[v]
        nodes.append(d)
    for e in gv.edges():
        source = idx[int(e.source())]
        target = idx[int(e.target())]
        strees = gv.edge_strees[e]
        d = dict(source=source, target=target, strees = list(strees),
                 taxedge=bool(gv.edge_in_taxonomy[e]))
        d['color'] = ecolor[e]
        d['width'] = ewidth[e]
        links.append(d)

    return dict(nodes=nodes, links=links)

def nexson2ptag(nexson):
    d = json.loads(nexson, encoding='utf-8')['nexml']
    otu_id2data = {}
    rv = {}
    for k,v in d['otusById'].iteritems():
        for otuid, otudata in v['otuById'].iteritems():
            assert otuid not in otu_id2data
            otu_id2data[otuid] = otudata

    treeids = []
    for k,v in d['treesById'].iteritems():
        for treeid, treedata in v['treeById'].iteritems():
            treeids.append(treeid)
            treedata['treeid'] = treeid
            r = buildtree(treedata, otu_id2data)
            lvs = r.leaves()
            prop = float(len([ lf for lf in lvs if lf.taxid ]))/len(lvs)
            if prop > 0.8 and len(lvs) < 5000:
                rv[treeid] = proctree(r)
    rv['treeids'] = treeids
    return rv
