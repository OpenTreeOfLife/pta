# -*- coding: utf-8 -*-
import os, treegraph

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
