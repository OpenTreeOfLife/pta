import sys, os, urllib2, subprocess, argparse, tarfile, tempfile, shutil, gzip
import codecs, sqlite3, json, ivy
import numpy as np
import ivy.treegraph as tg
gt = tg.gt # graph_tool.all
from glob import glob
from cStringIO import StringIO
from os.path import dirname, join, isdir
from git import Git
from itertools import chain

PTA_BASE = '/home/rree/src/pta'
PHYLESYSTEM_SHARDS = '/home/rree/src/phylesystem/shards'
OTT_GRAPH_FILENAME = 'ott2.8.gt.gz'
PTAG_DIR = join(PTA_BASE, 'static', 'ptag')

ap = argparse.ArgumentParser()
ap.add_argument('--fetch-ptags', dest='fetch_ptags', action='store_true')
ap.add_argument('--fetch-db', dest='fetch_db', action='store_true')
args = ap.parse_args()

cmd = 'find {} -name *.json'
nexsons = []
for d in [ join(PHYLESYSTEM_SHARDS, x)
           for x in os.listdir(PHYLESYSTEM_SHARDS)
           if isdir(join(PHYLESYSTEM_SHARDS, x)) ]:
    study_files = sorted(
        [ x for x in subprocess.check_output(cmd.format(d).split())
          .strip().split('\n') ])
    nexsons.append((Git(d), study_files))

def load_ott_graph():
    staticdir = join(PTA_BASE, 'static')
    localpath = join(staticdir, OTT_GRAPH_FILENAME)
    if not os.path.exists(localpath):
        remote = ('https://dl.dropboxusercontent.com/u/1939917/'+
                  OTT_GRAPH_FILENAME)
        cmd = 'curl -o {} {}'.format(localpath, remote)
        subprocess.check_call(cmd.split(), stdout=sys.stdout, stderr=sys.stderr)
    return tg.load_taxonomy_graph(localpath)

def proctree(t, otu_id2data):
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

def graph_json(g, dist=None, pos=None, ecolor=None, ewidth=None,
               vcolor=None, vsize=None, vtext=None, fp=None):
    nodes = []
    links = []
    idx = {}

    if pos:
        xmin = min([ pos[x][0] for x in g.vertices() ])
        ymin = min([ pos[x][1] for x in g.vertices() ])
        for x in g.vertices(): pos[x] = [pos[x][0]-xmin, pos[x][1]-ymin]

    for i,v in enumerate(g.vertices()):
        idx[int(v)] = i
        taxid = g.vertex_taxid[v]
        name = g.taxid_name(taxid) if taxid else ''#'node%s' % int(v)
        isleaf = v.out_degree()==0
        d = dict(label=name, isleaf=isleaf, strees=list(g.vertex_strees[v]))
        if taxid: d['taxid'] = taxid
        if dist: d['dist'] = dist[v]
        if pos and pos[v]:
            x, y = pos[v]
            d['x'] = x; d['y'] = y
            d['fixed'] = True
        if vcolor: d['color'] = vcolor[v]
        if vsize: d['size'] = vsize[v]
        if vtext: d['label'] = vtext[v]
        d['altlabel'] = g.vertex_name[v]
        nodes.append(d)
    for e in g.edges():
        source = idx[int(e.source())]
        target = idx[int(e.target())]
        strees = g.edge_strees[e]
        d = dict(source=source, target=target, strees = list(strees),
                 taxedge=bool(g.edge_in_taxonomy[e]))
        if ecolor: d['color'] = ecolor[e]
        if ewidth: d['width'] = ewidth[e]
        links.append(d)
    if fp:
        json.dump(dict(nodes=nodes, links=links), fp, indent=3)
    else:
        return json.dumps(dict(nodes=nodes, links=links), indent=3)

def make_graph_json(g, r, i, outfname):
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
                assert taxg.vertex_index[taxv]!=taxg.vertex_index[lf.v]
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

    with gzip.open(outfname, 'wb') as f:
        graph_json(gv, pos=pos, ecolor=ecolor, ewidth=ewidth,
                   vcolor=vcolor, vsize=vsize, vtext=gv.vertex_name, fp=f)

def index_taxa(fresh=False):
    g = load_ott_graph()
    dbfilename = 'phylesystem-taxon-index.db'
    dbdir = join(PTA_BASE, 'databases')
    dbpath = join(dbdir, dbfilename)
    if not os.path.exists(dbpath):
        fresh = True
    con = sqlite3.connect(dbpath)
    if fresh:
        cur = con.cursor()
        cur.execute('drop table if exists main')
        cur.execute(
            ('create table main ( '
             'id integer primary key, '
             'studyid text not null, '
             'mtime integer not null, '
             'citation text not null, '
             'treeid text not null, '
             'taxid integer not null, '
             'name text not null, '
             'monophyly_status integer not null default 0, '
             'otu integer not null default 0, '
             'tree_mrca integer not null default 0 '
             ')'))
    cur = con.cursor()
    fields = 'studyid, mtime, citation, treeid, taxid, name, otu, tree_mrca'
    INS = 'insert into main ({}) values (?,?,?,?,?,?,?,?)'.format(fields)
    for git, studies in nexsons:
        for x in studies:
            mtime = int(git.log('--format=%ct', x).split('\n')[0])
            path, fn = os.path.split(x)
            fnbase = os.path.splitext(fn)[0]
            studyid = fnbase
            if not fresh:
                sql = ('select count(*) from main where '
                       '(studyid = ?) and (mtime = ?)')
                if cur.execute(sql, (studyid, mtime)).fetchone()[0]:
                    print 'ignoring', fn
                    continue
            print 'indexing', fn#, mtime
            with codecs.open(x, encoding='utf-8') as f:
                try:
                    d = json.load(f)['nexml']
                except KeyError:
                    print fn, ': not a nexson file?'
                    continue
                try:
                    citation = d['^ot:studyPublicationReference']
                except KeyError:
                    print fn, ': no citation'
                    continue
                otu_id2data = {}
                for k,v in d['otusById'].iteritems():
                    for otuid, otudata in v['otuById'].iteritems():
                        assert otuid not in otu_id2data
                        otu_id2data[otuid] = otudata

                for k,v in d['treesById'].iteritems():
                    for treeid, treedata in v['treeById'].iteritems():
                        outfname = '{}.{}.{}.ptag.json.gz'.format(
                            fnbase, treeid, mtime)
                        p = join(PTAG_DIR, outfname)
                        if not os.path.exists(p):
                            continue
                        print 'inserting', treeid
                        nexecs = 0
                        treedata['nexson_file'] = x
                        treedata['treeid'] = treeid
                        r = proctree(treedata, otu_id2data)
                        lvs = r.leaves()
                        leaf_taxids = set([ x.taxid for x in lvs if x.taxid ])
                        if len(leaf_taxids) < 3:
                            print (treeid, ': only', len(leaf_taxids),
                                   'taxids found, skipping')
                            continue
                        for taxid in leaf_taxids:
                            name = g.taxid_unique_name(taxid)
                            row = (fnbase, #study_id
                                   mtime, citation, treeid, taxid, name, 1, 0)
                            try:
                                cur.execute(INS, row)
                                nexecs += 1
                            except sqlite3.IntegrityError:
                                print fn, ': can\'t execute', row
                                pass
                        rps = [ tg.taxid_rootpath(g, x) for x in leaf_taxids
                                if g.taxid_vertex.get(x) ]
                        mrca = tg.rootpath_mrca(rps)
                        name = g.taxid_unique_name(mrca)
                        row = (fnbase, mtime, citation, treeid, mrca, name, 0, 1)
                        try:
                            cur.execute(INS, row)
                            nexecs += 1
                        except sqlite3.IntegrityError:
                            print fn, ': can\'t execute', row
                            pass
                        internal_taxids = set()
                        for rp in rps:
                            for x in rp[1:]:
                                if x == mrca: break
                                internal_taxids.add(x)
                        for taxid in internal_taxids:
                            name = g.taxid_unique_name(taxid)
                            row = (fnbase, #study_id
                                   mtime, citation, treeid, taxid, name, 0, 0)
                            try:
                                cur.execute(INS, row)
                                nexecs += 1
                            except sqlite3.IntegrityError:
                                print fn, ': can\'t execute', row
                                pass
                print 'finished:', nexecs
                con.commit()
    if fresh:
        cur.execute('create index studyid_idx on main (studyid)')
        cur.execute('create index citation_idx on main (citation)')
        cur.execute('create index treeid_idx on main (treeid)')
        cur.execute('create index taxid_idx on main (taxid)')
        cur.execute('create index name_idx on main (name)')
        cur.execute('create index mono_idx on main (monophyly_status)')
        cur.execute('create index otu_idx on main (otu)')
        cur.execute('create index tree_mrca_idx on main (tree_mrca)')

def make_ptags():
    g = load_ott_graph()
    for git, studies in nexsons:
        for x in studies:
            mtime = git.log('--format=%ct', x).split('\n')[0]
            path, fn = os.path.split(x)
            fnbase = os.path.splitext(fn)[0]
            with codecs.open(x, encoding='utf-8') as f:
                try:
                    d = json.load(f)['nexml']
                except KeyError:
                    print fn, ': not a nexson file?'
                    continue
                otu_id2data = {}
                for k,v in d['otusById'].iteritems():
                    for otuid, otudata in v['otuById'].iteritems():
                        assert otuid not in otu_id2data
                        otu_id2data[otuid] = otudata

                for k,v in d['treesById'].iteritems():
                    for treeid, treedata in v['treeById'].iteritems():
                        outfname = '{}.{}.{}.ptag.json.gz'.format(
                            fnbase, treeid, mtime)
                        p = join(PTAG_DIR, outfname)
                        if not os.path.exists(p):
                            treedata['nexson_file'] = x
                            treedata['treeid'] = treeid
                            try:
                                r = proctree(treedata, otu_id2data)
                                assert r
                            except:
                                print 'ROOT FAIL:', fn, treeid
                                continue
                            lvs = r.leaves()
                            prop = float(len([ lf for lf in lvs
                                               if lf.taxid ]))/len(lvs)
                            if prop > 0.8 and len(lvs) < 5000:
                                try:
                                    make_graph_json(g, r, 1, p)
                                    print outfname
                                except:
                                    print 'FAIL:', fn, treeid
                                    

def taxon_freqs():
    from collections import Counter
    c = Counter()
    for fname in glob(join(PTAG_DIR, '*.ptag.json')):
        with open(fname) as f:
            d = json.load(f)
        for n in d['nodes']:
            t = n.get('label'); tid = n.get('taxid')
            if (not tid) or (tid and t and len(t.split()))>1:
                continue
            c[(t,tid)] += 1
    with open(join(PTA_BASE, 'static', 'taxon-freqs.json'), 'w') as f:
        json.dump([ dict(name=k[0],freq=v,taxid=k[1]) for k,v in
                    c.most_common(100) ], f, indent=2)

def main():
    dbfilename = 'phylesystem-taxon-index.db'
    dbdir = join(PTA_BASE, 'databases')
    dbpath = join(dbdir, dbfilename)
    staticdir = join(PTA_BASE, 'static')
    ptagdir = join(staticdir, 'ptag')
    
    remotedb = 'https://dl.dropboxusercontent.com/u/1939917/'+dbfilename

    if (not os.path.isfile(dbpath)) or args.fetch_db:
        print 'fetching database...'
        sys.stdout.flush()
        cmd = 'curl -o {} {}'.format(dbpath, remotedb)
        subprocess.check_call(cmd.split(), stdout=sys.stdout, stderr=sys.stderr)

    if not os.path.exists(ptagdir):
        print 'ptag directory not found; creating {}'.format(ptagdir)
        os.mkdir(ptagdir)

    local_ptags = glob(join(ptagdir,'*.ptag.json'))
    
    if (not local_ptags) or args.fetch_ptags:
        print 'fetching remote ptags...',
        remote = 'https://dl.dropboxusercontent.com/u/1939917/ptag.tar.gz'
        sys.stdout.flush()
        r = StringIO(urllib2.urlopen(remote).read())
        tf = tarfile.open(fileobj=r, mode='r:gz')
        d = tempfile.mkdtemp()
        tf.extractall(d)
        cmd = 'rsync -a {} {}'.format(join(d,'ptag'), staticdir)
        subprocess.check_call(cmd.split(), stdout=sys.stdout, stderr=sys.stderr)
        shutil.rmtree(d)

def mismapped_leaves(r):
    # r is a mapped tree from tg.map_stree
    lvs = r.leaves()
    thresh = 3.5
    maxdepth = float(max([ len(list(lf.rootpath())) for lf in lvs ]))
    a = np.array([ len(lf.taxids)/maxdepth for lf in lvs ])
    v = np.subtract(a, np.median(a)).clip(min=0)
    med = np.median(v[v>0])
    return [ lvs[i] for i in np.nonzero((0.6745 * v/med) > thresh)[0] ]
        
def index_tree(g, cur, r, studyid, mtime, citation, treeid):
    # r is the root node after proctree
    lvs = r.leaves()
    rps = [ tg.taxid_rootpath(g, lf.taxid) for lf in lvs if lf.taxid ]
    mrca = tg.rootpath_mrca(rps)
    rps = [ p[:p.index(mrca)+1] for p in rps ]
    taxids = set(chain.from_iterable(rps))

    taxg = tg.taxid_new_subgraph(g, taxids)
    tg.map_stree(taxg, r)
    verts = taxg.new_vertex_property('bool')
    edges = taxg.new_edge_property('bool')
    # add stree's nodes and branches into taxonomy graph
    tg.merge_stree(taxg, r, 1, verts, edges)
    v2node = dict([ (n.v, n) for n in r ])

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
                assert taxg.vertex_index[taxv]!=taxg.vertex_index[lf.v]
                e = taxg.add_edge(taxv, lf.v)
            taxg.edge_in_taxonomy[e] = 1

    fields = ('studyid, mtime, citation, treeid, taxid, name, '
              'monophyly_status, otu, tree_mrca')
    INS = 'insert into main ({}) values (?,?,?,?,?,?,?,?,?)'.format(fields)

    # assess monophyly
    vname = taxg.vp['name']
    uname = taxg.vp['unique_name']
    for v in taxg.vertices():
        taxid = taxg.vertex_taxid[v]
        if taxid:
            name = uname[v] or vname[v]
            traced = bool(taxg.vertex_strees[v])
            odeg = v.out_degree()

            mono = 0; otu = 0; mrca = 0

            if not traced:
                mono = -1
            if traced and odeg > 1:
                mono = 1

            if v.in_degree()==0: # root
                mrca = 1

            if v.out_degree()==0:
                otu = 1
        
            x = (studyid, mtime, citation, treeid, taxid, name, mono, otu, mrca)
            try:
                cur.execute(INS, x)
            except sqlite3.IntegrityError:
                print studyid, treeid, ': can\'t execute', x

def experiment():
    g = load_ott_graph()
    ## x = join(PHYLESYSTEM_SHARDS, 'phylesystem-1/study/pg_14/pg_14/pg_14.json')
    ## treeid = 'tree12'
    x = join(PHYLESYSTEM_SHARDS,
             'phylesystem-1/study/pg_71/pg_2071/pg_2071.json')
    treeid = 'tree4264'
    with codecs.open(x, encoding='utf-8') as f:
        d = json.load(f)['nexml']
    otu_id2data = {}
    for k,v in d['otusById'].iteritems():
        for otuid, otudata in v['otuById'].iteritems():
            assert otuid not in otu_id2data
            otu_id2data[otuid] = otudata

    for dv in d['treesById'].values():
        if treeid in dv['treeById']:
            break

    treedata = dv['treeById'][treeid]
    treedata['treeid'] = treeid
    r = proctree(treedata, otu_id2data)
    lvs = r.leaves()
            
    rps = [ tg.taxid_rootpath(g, lf.taxid) for lf in lvs if lf.taxid ]
    mrca = tg.rootpath_mrca(rps)
    rps = [ p[:p.index(mrca)+1] for p in rps ]
    taxids = set(chain.from_iterable(rps))

    taxg = tg.taxid_new_subgraph(g, taxids)
    tg.map_stree(taxg, r)
    verts = taxg.new_vertex_property('bool')
    edges = taxg.new_edge_property('bool')

    # add stree's nodes and branches into taxonomy graph
    tg.merge_stree(taxg, r, 1, verts, edges)

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
                assert taxg.vertex_index[taxv]!=taxg.vertex_index[lf.v]
                e = taxg.add_edge(taxv, lf.v)
            taxg.edge_in_taxonomy[e] = 1

    # assess monophyly
    vname = taxg.vp['name']
    uname = taxg.vp['unique_name']
    for v in taxg.vertices():
        if v.in_degree()==0:
            continue
        tid = taxg.vertex_taxid[v]
        if tid:
            print v, tid, uname[v] or vname[v], bool(taxg.vertex_strees[v]), v.out_degree()
