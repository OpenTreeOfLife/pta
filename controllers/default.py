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
    
    if not d:
        d = request.vars.data
    if d:
        data = treegraph.nexson2ptag(d)
    if data:
        print data.get('treeids')
    return dict(data=data, form=form)
