# -*- coding: utf-8 -*-
import os, treegraph, requests, json
from gzip import GzipFile
from collections import defaultdict

def index():
    d = data = None
    form = SQLFORM.factory(
        Field('data', 'upload', uploadfolder='/tmp')
        )
    if form.process().accepted:
        fp = os.path.join('/tmp', form.vars.data)
        d = open(fp).read()

    try:
        session.filename = request.vars.data.filename
    except:
        session.filename = None
    
    if not d:
        d = request.vars.data
    if d:
        data = treegraph.nexson2ptag(d)
    if data:
        for treeid in data.get('treeids'):
            k = 'pta_{}'.format(treeid)
            session[k] = data[treeid]
    return dict(data=data, form=form)

def view():
    treeid = request.vars.treeid
    k = 'pta_{}'.format(treeid)
    d = session.get(k)
    return dict(data=d)

def _ptag_file_exists(studyid, treeid, mtime):
    ptag = '{}/static/ptag/{}.{}.{}.ptag.json.gz'.format(
        request.folder, studyid, treeid, mtime)
    return os.path.exists(ptag)

def view2():
    studyid, treeid = request.args
    ts = treegraph.study_timestamp(studyid)
    ptag = '{}/static/ptag/{}.{}.{}.ptag.json.gz'.format(
        request.folder, studyid, treeid, ts)
    if os.path.exists(ptag):
        data = XML(GzipFile(ptag).read())
    else:
        u = 'http://api.opentreeoflife.org/v2/study/{}.json'
        p = dict(output_nexml2json='1.2.1')
        r = requests.get(u.format(studyid), params=p)
        data = XML(json.dumps(treegraph.nexson2ptag(r.text, treeid)[treeid]))
    return dict(data=data)

def search():
    g = treegraph.g
    search_option = request.vars.search_option or '0'
    search_term = request.vars.search_term or ''
    t = db.main
    f = t.name
    rv = []; cit2trees = defaultdict(list)

    if len(request.args)==2 and request.args[0]=='taxid':
        tid = int(request.args[1])
        name = g.taxid_unique_name(tid)
        search_term = name
        rows = db(f==name).select(t.studyid, t.treeid, t.citation, t.mtime,
                                  distinct=True)
        for r in rows:
            if _ptag_file_exists(r.studyid, r.treeid, r.mtime):
                q = ((t.studyid==r.studyid)&
                     (t.treeid==r.treeid)&
                     (t.tree_mrca==1))
                n = db(q).select(t.name, limitby=(0,1)).first().name
                rv.append((r, n))
        for r, n in rv:
            cit2trees[r.citation].append((r.treeid, n, r.studyid))

    form = SQLFORM.factory(
        Field('search_option', requires=IS_IN_SET((
            ('0', 'represented anywhere in the tree'),
            ('1', 'an OTU mapped to a leaf node(s)'),
            ('2', 'the MRCA of all leaves')
            ), zero=None), label='where taxon is',
            default=search_option),
        Field('search_term', requires=IS_NOT_EMPTY(), label='and named',
              default=search_term),
        submit_button='Find trees',
        table_name='myform'
        )

    if form.process(message_onsuccess=None, keepvalues=True).accepted:
        term = form.vars.search_term
        opt = form.vars.search_option or '0'
        if term:
            print term
            q = f.like(term)
            if opt == '1':
                q &= (t.otu==1)
            elif opt == '2':
                q &= (t.tree_mrca==1)
            rows = db(q).select(t.studyid, t.treeid, t.citation, t.mtime,
                                distinct=True)
            for r in rows:
                if _ptag_file_exists(r.studyid, r.treeid, r.mtime):
                    q = ((t.studyid==r.studyid)&
                         (t.treeid==r.treeid)&
                         (t.tree_mrca==1))
                    n = db(q).select(t.name, limitby=(0,1)).first().name
                    rv.append((r, n))
            for r, n in rv:
                cit2trees[r.citation].append((r.treeid, n, r.studyid))
    return dict(form=form, rows=rv, cit2trees=cit2trees)

def taxon_cloud():
    response.files.append("http://d3js.org/d3.v3.min.js")
    response.files.append(URL('static', 'js/d3.layout.cloud.js'))
    d = json.load(open('{}/static/taxon-freqs.json'.format(request.folder)))
    j = json.dumps([ dict(name=x['name'], freq=x['freq'], taxid=x['taxid'],
                          url=URL('search', args=['taxid', x['taxid']],
                                  vars=dict(search_option=0)))
                     for x in d ])
    return dict(data=j)

def name_search_autocomplete():
    rv = []
    term = request.vars.term or ''
    opt = request.vars.search_option or '0'
    t = db.main
    f = t.name
    if term:
        q = f.like('{}%'.format(term))
        if opt == '1':
            q &= (t.otu==1)
        elif opt == '2':
            q &= (t.tree_mrca==1)
        rows = db(q).select(f, distinct=True, orderby=f, limitby=(0,25))
        rv = [ x.name for x in rows ]
    return response.json(rv)



