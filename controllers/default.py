# -*- coding: utf-8 -*-
import os, treegraph, requests, json
from gzip import GzipFile

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
    form = SQLFORM.factory(Field('term', 'string'))
    t = db.main
    f = t.name
    rv = []
    if form.process(message_onsuccess=None).accepted:
        s = form.vars.term
        if s:
            q = f.like(s) & (t.otu==1)
            rows = db(q).select(t.studyid, t.treeid, t.citation, t.mtime,
                                distinct=True)
            for r in rows:
                if _ptag_file_exists(r.studyid, r.treeid, r.mtime):
                    q = ((t.studyid==r.studyid)&
                         (t.treeid==r.treeid)&
                         (t.tree_mrca==1))
                    n = db(q).select(t.name, limitby=(0,1)).first().name
                    rv.append((r, n))
    return dict(form=form, rows=rv)

def name_search_autocomplete():
    rv = []
    term = request.vars.term or ''
    t = db.main
    f = t.name
    if len(term) > 2:
        q = f.like('%{}%'.format(term)) & (t.otu==1)
        rows = db(q).select(f, distinct=True, orderby=f, limitby=(0,25))
        rv = [ x.name for x in rows ]
    return response.json(rv)
