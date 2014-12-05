import ivy, json
from peyotl.api.phylesystem_api import PhylesystemAPI
psapi = PhylesystemAPI(get_from='local').phylesystem_obj
import git
from ivy import treegraph as tg
gt = tg.gt

g = tg.load_taxonomy_graph('applications/pta/static/ott2.8.gt.gz')

def buildtree(t, otu_id2data):
    from ivy.tree import Node
    if t.has_key('nodeById'):
        # newer Nexson 1.2.1
        node_id2data = t['nodeById']
    else:
        # older Nexson 1.0.0
        node_id2data = {} 
        for n in t['node']:
            node_id2data[ n['@id'] ] = n
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
    if t.has_key('edgeBySourceId'):
        # newer Nexson 1.2.1
        edges_by_source_id = t['edgeBySourceId']
    else:
        # older Nexson 1.0.0
        edges_by_source_id = {} 
        for e in t['edge']:
            if not edges_by_source_id.has_key( e['@source'] ):
                edges_by_source_id[ e['@source'] ] = {}
            edges_by_source_id[ e['@source'] ][ e['@id'] ] = e
    for nid, ed in edges_by_source_id.iteritems():
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
        try:
            name = gv.vertex_name[v]
        except:
            name = gv.taxid_name(taxid) if taxid else ''#'node%s' % int(v)
        isleaf = v.out_degree()==0
        d = dict(label=name, isleaf=isleaf, strees=list(gv.vertex_strees[v]),
                 altlabel=name)
        if taxid: d['taxid'] = taxid
        ## if dist: d['dist'] = dist[v]
        if pos and pos[v]:
            x, y = pos[v]
            d['x'] = x; d['y'] = y
            d['fixed'] = True
        d['color'] = vcolor[v]
        d['size'] = vsize[v]
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

def nexson2ptag(nexson, only_treeid=None):
    # nexson is a string
    try:
        d = json.loads(nexson, encoding='utf-8')['nexml']
    except:
        d = json.loads(nexson, encoding='utf-8')['data']['nexml']
    otu_id2data = {}
    rv = {}
    treeids = []
    if d.has_key('otusById'):
        # we're using a newer version of Nexson, probably 1.2.1
        for k,v in d['otusById'].iteritems():
            for otuid, otudata in v['otuById'].iteritems():
                assert otuid not in otu_id2data
                otu_id2data[otuid] = otudata
        for k,v in d['treesById'].iteritems():
            for treeid, treedata in v['treeById'].iteritems():
                if only_treeid and treeid != treeid:
                    continue
                treeids.append(treeid)
                treedata['treeid'] = treeid
                r = buildtree(treedata, otu_id2data)
                lvs = r.leaves()
                prop = float(len([ lf for lf in lvs if lf.taxid ]))/len(lvs)
                if prop > 0.8 and len(lvs) < 5000:
                    rv[treeid] = proctree(r)
    else:
        # this is older Nexson (1.0.0), possibly from the curation app
        for otus_collection in d['otus']:
            print otus_collection
            for otudata in otus_collection['otu']:
                otuid = otudata['@id']
                assert otuid not in otu_id2data
                otu_id2data[otuid] = otudata
        for trees_collection in d['trees']:
            for treedata in trees_collection['tree']:
                treeid = treedata['@id']
                if only_treeid and treeid != treeid:
                    continue
                treeids.append(treeid)
                treedata['treeid'] = treeid
                r = buildtree(treedata, otu_id2data)
                lvs = r.leaves()
                prop = float(len([ lf for lf in lvs if lf.taxid ]))/len(lvs)
                if prop > 0.8 and len(lvs) < 5000:
                    rv[treeid] = proctree(r)

    rv['treeids'] = treeids
    return rv

def study_timestamp(studyid):
    repo = psapi.get_shard(studyid).path
    g = git.Git(repo)
    pth = psapi.get_filepath_for_study(studyid)
    return git.log('--format=%ct', pth).split('\n')[0]
