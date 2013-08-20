#!/usr/bin/env python3
# coding: utf-8
#
# change_log.py

import parsers
import urllib
import os
import write_log
import datetime
import pickle
import ipdb
import parsed
import time
from changelog_config import changelog_data_opt
from constants import RES_LOCATION, PARSER_TYPE

# parse any files needed for change-log resolution
# for label, data_tuple in changelog_data_opt.items():
#     url = data_tuple[RES_LOCATION]
#     parser = data_tuple[PARSER_TYPE](url)
#     print('Running ' +str(parser))
#     for x in parser.parse():
#         parsed.build_data(x, str(parser))


# BELNamespaceParser - parse() returns the specific url for each namespace
# currently published on resource.belframework.org
start_time = time.time()
parser = parsers.BELNamespaceParser()
print('Running BELNamespace_Parser')
old_entrez = set()
old_hgnc = set()
old_mgi = set()
old_rgd = set()
old_sp_names = set()
old_sp_ids = set()
old_affy = set()
old_chebi_names = set()
old_chebi_ids = set()
old_gobp_ns_ids = set()
old_gobp_ns_names = set()
old_gocc_ns_ids = set()
old_gocc_ns_names = set()
old_mesh_bio = set()
old_mesh_cell = set()
old_mesh_diseases = set()

old_chebi_eq_names = dict()
old_chebi_eq_ids = dict()
old_gobp_eq_ids = dict()
old_gobp_eq_names = dict()
old_gocc_eq_ids = dict()
old_gocc_eq_names = dict()

# iterate over the urls to the .belns files, collecting the entries
# from the old data.
for url in parser.parse():
    debug = False
    namespaces = { 'entrez' : (False, old_entrez),
                   'hgnc' : (False, old_hgnc),
                   'mgi' : (False, old_mgi),
                   'rgd' : (False, old_rgd),
                   'swissprot-entry' : (False, old_sp_names),
                   'swissprot-accession' : (False, old_sp_ids),
                   'affy' : (False, old_affy),
                   'chebi-name' : (False, old_chebi_names),
                   'chebi-id' : (False, old_chebi_ids),
                   'go-biological-processes-acc' : (False, old_gobp_ns_ids),
                   'go-biological-processes-names' : (False, old_gobp_ns_names),
                   'go-cellular-component-terms' : (False, old_gocc_ns_names),
                   'go-cellular-component-acc' : (False, old_gocc_ns_ids),
                   'mesh-bio' : (False, old_mesh_bio),
                   'mesh-cell' : (False, old_mesh_cell),
                   'mesh-diseases' : (False, old_mesh_diseases)}

    open_url = urllib.request.urlopen(url)
    for ns in namespaces:
        if ns in open_url.url:
            namespaces[ns] = (True, namespaces[ns][1])
    marker = False
    for u in open_url:
        if '[Values]' in str(u):
            marker = True
            continue
        if marker is False:
            continue
        # we are into namespace pairs with '|' delimiter
        t = u.decode('utf-8')
        tokenized = t.split('|')
        token = tokenized[0]
        for k, v in namespaces.items():
            if v[0]:
                v[1].add(token)

# parse the old .beleq files needed for resolving lost values
parser = parsers.BELEquivalenceParser()
print('Running BELEquivalence_Parser')
for url in parser.parse():

    equivalences = {
        'go-biological-processes-acc' : (False, old_gobp_eq_ids),
        'go-biological-processes-names' : (False, old_gobp_eq_names),
        'go-cellular-component-terms' : (False, old_gocc_eq_names),
        'go-cellular-component-acc' : (False, old_gocc_eq_ids),
        'chebi-names' : (False, old_chebi_eq_names),
        'chebi-ids' : (False, old_chebi_eq_ids)}

    open_url = urllib.request.urlopen(url)
    for eq in equivalences:
        if eq in open_url.url:
            equivalences[ns] = (True, equivalences[eq][1])
    marker = False
    for u in open_url:
        if '[Values]' in str(u):
            marker = True
            continue
        if marker is False:
            continue
        # we are into namespace pairs with '|' delimiter
        t = u.decode('utf-8')
        tokenized = t.split('|')
        value = tokenized[0]
        uid = tokenized[1]
        for k, v in equivalences.items():
            if v[0]:
                v[1][value] = uid

# reverse the dicts, the uuids will be the keys
inv_gobp_ids = {v:k for k, v in old_gobp_eq_ids.items()}
inv_gobp_names = {v:k for k, v in old_gobp_eq_names.items()}
inv_gocc_ids = {v:k for k, v in old_gocc_eq_ids.items()}
inv_gocc_names = {v:k for k, v in old_gocc_eq_names.items()}

inv_chebi_ids = {v:k for k, v in old_chebi_eq_ids.items()}
inv_chebi_names = {v:k for k, v in old_chebi_eq_names.items()}

gobp_names_to_ids = {}
gocc_names_to_ids = {}
chebi_names_to_ids = {}
# any keys that match up, get the corresponding name/id
for uid in inv_gobp_ids:
    if uid in inv_gobp_names:
        name = inv_gobp_names.get(uid)
        id = inv_gobp_ids.get(uid)
        gobp_names_to_ids[name] = id
for uid in inv_gocc_ids:
    if uid in inv_gocc_names:
        name = inv_gocc_names.get(uid)
        id = inv_gocc_ids.get(uid)
        gocc_names_to_ids[name] = id
for uid in inv_chebi_ids:
    if uid in inv_chebi_names:
        name = inv_chebi_names.get(uid)
        id = inv_chebi_ids.get(uid)
        chebi_names_to_ids[name] = id

def name_to_id(name):
    if name in gobp_names_to_ids:
        return gobp_names_to_ids[name]
    elif name in gocc_names_to_ids:
        return gocc_names_to_ids[name]
    elif name in chebi_names_to_ids:
        return chebi_names_to_ids[name]
    else:
        return None

print('===========================================')
print('len of old entrez is ' +str(len(old_entrez)))
print('len of old hgnc is ' +str(len(old_hgnc)))
print('len of old mgi is ' +str(len(old_mgi)))
print('len of old rgd is ' +str(len(old_rgd)))
print('len of old swissprot names is ' +str(len(old_sp_names)))
print('len of old swissprot accessions is ' +str(len(old_sp_ids)))
print('len of old affy is ' +str(len(old_affy)))
print('len of old chebi names is ' +str(len(old_chebi_names)))
print('len of old chebi ids is ' +str(len(old_chebi_ids)))
print('len of old gobp names is ' +str(len(old_gobp_ns_names)))
print('len of old gobp ids is ' +str(len(old_gobp_ns_ids)))
print('len of old gocc names is ' +str(len(old_gocc_ns_names)))
print('len of old gocc ids is ' +str(len(old_gocc_ns_ids)))
print('len of old mesh-bio is ' +str(len(old_mesh_bio)))
print('len of old mesh-cell is ' +str(len(old_mesh_cell)))
print('len of old mesh-diseases is ' +str(len(old_mesh_diseases)))
print('===========================================')

new_entrez = set()
new_hgnc = set()
new_mgi = set()
new_rgd = set()
new_sp_names = set()
new_sp_ids = set()
new_affy = set()
new_chebi_names = set()
new_chebi_ids = set()
new_gobp_ns_ids = set()
new_gobp_ns_names = set()
new_gocc_ns_ids = set()
new_gocc_ns_names = set()
new_mesh_bio = set()
new_mesh_cell = set()
new_mesh_diseases = set()

# gather the new data for comparison (locally stored for now)
indir = '/home/jhourani/openbel-contributions/resource_generator/lions'
for root, dirs, filenames in os.walk(indir):
    for f in filenames:
        if '.belns' in f:
            with open(os.path.join(root, f), 'r') as fp:
                if 'entrez' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_entrez.add(token)
                elif 'hgnc' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_hgnc.add(token)
                elif 'mgi' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_mgi.add(token)
                elif 'rgd' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_rgd.add(token)
                elif 'swiss-acc' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_sp_ids.add(token)
                elif 'swiss-entry' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_sp_names.add(token)
                elif 'affy' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_affy.add(token)
                elif 'chebi-names' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_chebi_names.add(token)
                elif 'chebi-ids' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_chebi_ids.add(token)
                elif 'go-biological-processes-acc' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_gobp_ns_ids.add(token)
                elif 'go-biological-processes-names' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_gobp_ns_names.add(token)
                elif 'go-cellular-component-acc' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_gocc_ns_ids.add(token)
                elif 'go-cellular-component-terms' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_gocc_ns_names.add(token)
                elif 'mesh-bio' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_mesh_bio.add(token)
                elif 'mesh-cell' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_mesh_cell.add(token)
                elif 'mesh-diseases' in fp.name:
                    for line in fp:
                        # if '[Values]' in str(line):
                        #     marker = True
                        #     continue
                        # if marker is False:
                        #     continue
                        tokenized = str(line).split('|')
                        token = tokenized[0]
                        new_mesh_diseases.add(token)

print('len of new entrez is ' +str(len(new_entrez)))
print('len of new hgnc is ' +str(len(new_hgnc)))
print('len of new mgi is ' +str(len(new_mgi)))
print('len of new rgd is ' +str(len(new_rgd)))
print('len of new swiss names is ' +str(len(new_sp_names)))
print('len of new swiss ids is ' +str(len(new_sp_ids)))
print('len of new affy is ' +str(len(new_affy)))
print('len of new chebi names is ' +str(len(new_chebi_names)))
print('len of new chebi ids is ' +str(len(new_chebi_ids)))
print('len of new gobp names is ' +str(len(new_gobp_ns_names)))
print('len of new gobp ids is ' +str(len(new_gobp_ns_ids)))
print('len of new gocc names is ' +str(len(new_gocc_ns_names)))
print('len of new gocc ids is ' +str(len(new_gocc_ns_ids)))
print('len of new mesh-bio is ' +str(len(new_mesh_bio)))
print('len of new mesh-cell is ' +str(len(new_mesh_cell)))
print('len of new mesh-diseases is ' +str(len(new_mesh_diseases)))

# values in the old data that are not in the new (either withdrawn or replaced)
entrez_lost = [x for x in old_entrez if x not in new_entrez]
hgnc_lost = [x for x in old_hgnc if x not in new_hgnc]
mgi_lost = [x for x in old_mgi if x not in new_mgi]
rgd_lost = [x for x in old_rgd if x not in new_rgd]
sp_lost_names = [x for x in old_sp_names if x not in new_sp_names]
sp_lost_ids = [x for x in old_sp_ids if x not in new_sp_ids]
affy_lost = [x for x in old_affy if x not in new_affy]
chebi_lost_names = [x for x in old_chebi_names if x not in new_chebi_names]
chebi_lost_ids = [x for x in old_chebi_ids if x not in new_chebi_ids]
gobp_lost_ns_names = [x for x in old_gobp_ns_names if x not in new_gobp_ns_names]
gobp_lost_ns_ids = [x for x in old_gobp_ns_ids if x not in new_gobp_ns_ids]
gocc_lost_ns_names = [x for x in old_gocc_ns_names if x not in new_gocc_ns_names]
gocc_lost_ns_ids = [x for x in old_gocc_ns_ids if x not in new_gocc_ns_ids]
mesh_bio_lost = [x for x in old_mesh_bio if x not in new_mesh_bio]
mesh_cell_lost = [x for x in old_mesh_cell if x not in new_mesh_cell]
mesh_diseases_lost = [x for x in old_mesh_diseases if x not in new_mesh_diseases]

print('===========================================')
print('lost entrez values ' +str(len(entrez_lost)))
print('lost hgnc values ' +str(len(hgnc_lost)))
print('lost mgi values ' +str(len(mgi_lost)))
print('lost rgd values ' +str(len(rgd_lost)))
print('lost swissprot names ' +str(len(sp_lost_names)))
print('lost swissprot ids ' +str(len(sp_lost_ids)))
print('lost affy values ' +str(len(affy_lost)))
print('lost chebi names ' +str(len(chebi_lost_names)))
print('lost chebi ids ' +str(len(chebi_lost_ids)))
print('lost gobp names ' +str(len(gobp_lost_ns_names)))
print('lost gobp ids ' +str(len(gobp_lost_ns_ids)))
print('lost gocc names ' +str(len(gocc_lost_ns_names)))
print('lost gocc ids ' +str(len(gocc_lost_ns_ids)))
print('lost mesh_bio ' +str(len(mesh_bio_lost)))
print('lost mesh_cell ' +str(len(mesh_cell_lost)))
print('lost mesh_diseases ' +str(len(mesh_diseases_lost)))

# values in the new data that are not in the old (either new or a replacement)
entrez_gained = [x for x in new_entrez if x not in old_entrez]
hgnc_gained = [x for x in new_hgnc if x not in old_hgnc]
mgi_gained = [x for x in new_mgi if x not in old_mgi]
rgd_gained = [x for x in new_rgd if x not in old_rgd]
sp_gained_names = [x for x in new_sp_names if x not in old_sp_names]
sp_gained_ids = [x for x in new_sp_ids if x not in old_sp_ids]
affy_gained = [x for x in new_affy if x not in old_affy]
chebi_gained_names = [x for x in new_chebi_names if x not in old_chebi_names]
chebi_gained_ids = [x for x in new_chebi_ids if x not in old_chebi_ids]
gobp_gained_ns_names = [x for x in new_gobp_ns_names if x not in old_gobp_ns_names]
gobp_gained_ns_ids = [x for x in new_gobp_ns_ids if x not in old_gobp_ns_ids]
gocc_gained_ns_names = [x for x in new_gocc_ns_names if x not in old_gocc_ns_names]
gocc_gained_ns_ids = [x for x in new_gocc_ns_ids if x not in old_gocc_ns_ids]
mesh_gained_bio = [x for x in new_mesh_bio if x not in old_mesh_bio]
mesh_gained_cell = [x for x in new_mesh_cell if x not in old_mesh_cell]
mesh_gained_diseases = [x for x in new_mesh_diseases if x not in old_mesh_diseases]

print('===========================================')
print('gained entrez values ' +str(len(entrez_gained)))
print('gained hgnc values ' +str(len(hgnc_gained)))
print('gained mgi values ' +str(len(mgi_gained)))
print('gained rgd values ' +str(len(rgd_gained)))
print('gained swissprot names ' +str(len(sp_gained_names)))
print('gained swissprot ids ' +str(len(sp_gained_ids)))
print('gained affy values ' +str(len(affy_gained)))
print('gained chebi names ' +str(len(chebi_gained_names)))
print('gained chebi ids ' +str(len(chebi_gained_ids)))
print('gained gobp names ' +str(len(gobp_gained_ns_names)))
print('gained gobp ids ' +str(len(gobp_gained_ns_ids)))
print('gained gocc names ' +str(len(gocc_gained_ns_names)))
print('gained gocc ids ' +str(len(gocc_gained_ns_ids)))
print('gained mesh_bio ' +str(len(mesh_gained_bio)))
print('gained mesh_cell ' +str(len(mesh_gained_cell)))
print('gained mesh_diseases ' +str(len(mesh_gained_diseases)))
print('===========================================')

# with open('hgnc-new-values.txt', 'w') as fp:
#     for val in hgnc_gained:
#         fp.write(val +'\n')

# iterate lost values for each dataset, find out if they are withdrawn or
# replaced. If replaced, map oldname->newname. Otherwise map oldname->withdrawn.
change_log = {}
change_log['entrez'] = {}
change_log['hgnc'] = {}
change_log['mgi'] = {}
change_log['rgd'] = {}
change_log['swiss-names'] = {}
change_log['swiss-ids'] = {}
change_log['gobp-ids'] = {}
change_log['gobp-names'] = {}
change_log['gocc-ids'] = {}
change_log['gocc-names'] = {}
change_log['chebi-names'] = {}
change_log['chebi-ids'] = {}
change_log['mesh-bio'] = {}
change_log['mesh-cell'] = {}
change_log['mesh-diseases'] = {}
previous_symbols = []
previous_names = []
symbols_and_names = []
sp_accession_ids = []
for label, data_tuple in changelog_data_opt.items():
    url = data_tuple[RES_LOCATION]
    parser = data_tuple[PARSER_TYPE](url)

    if str(parser) == 'EntrezGeneHistory_Parser':
        print('Gathering Entrez update info...')
        with open('resolved-entrez.txt', 'w') as fp:
            for row in parser.parse():
                discontinued_id  = row.get('Discontinued_GeneID')
                fp.write(str(discontinued_id)+'\n')
                gid = row.get('GeneID')
                replacement_id = gid if gid != '-' else 'withdrawn'
                log = change_log.get('entrez')
                log[discontinued_id] = replacement_id

    elif str(parser) == 'HGNC_Parser':
        print('Gathering HGNC update info...')
        with open('unresolved-hgnc.txt', 'w') as fp:
            for row in parser.parse():
                val = row.get('Approved Symbol')
                # ps = row.get('Previous Symbols')
                # pn = row.get('Previous Names')
                # if len(ps) > 0:
                #     previous_symbols.extend(ps.split(', '))
                # if len(pn) > 0:
                #     previous_names.extend(pn.split(', '))
                if '~withdrawn' in val:
                    approved_name = row.get('Approved Name')
                    # no replacement
                    if 'entry withdrawn' in approved_name:
                        old_val = val.split('~')[0]
                        log = change_log.get('hgnc')
                        log[old_val] = 'withdrawn'
                    # has a replacement
                    if 'symbol withdrawn' in approved_name:
                        old_val = val.split('~')[0]
                        new_val = approved_name.split('see ')[1]
                        log = change_log.get('hgnc')
                        log[old_val] = new_val
            # symbols_and_names = previous_symbols + previous_names
            # unresolved = [x for x in hgnc_lost if x not in symbols_and_names \
            #                   and x not in change_log.get('hgnc')]
            # for u in unresolved:
            #     fp.write(str(u) +'\n')
            # print('Really unresolved: ' +str(len(unresolved)))

    elif str(parser) == 'MGI_Parser':
        print('Gathering MGI update info...')
        with open('resolved-mgi.txt', 'w') as fp:
            for row in parser.parse():
                mgi_accession = row.get('MGI Accession ID')
                if mgi_accession == 'NULL':
                    old_val = row.get('Marker Symbol')
                    name = row.get('Marker Name')
                    fp.write(str(old_val)+'\n')
                    if '=' in name:
                        log = change_log.get('mgi')
                        log[old_val] = name.split('= ')[1]
                    elif 'withdrawn' in name:
                        log = change_log.get('mgi')
                        log[old_val] = 'withdrawn'

    elif str(parser) == 'RGD_Parser':
        # rgd changes (still dont know if withdrawn or replaced!!)
        print('Gathering RGD update info...')
        with open('resolved-rgd.txt', 'w') as fp:
            for row in parser.parse():
                new_val = row.get('SYMBOL')
                lost_vals = row.get('OLD_SYMBOL').split(';')
                if len(lost_vals) != 0:
                    for symbol in lost_vals:
                        fp.write(str(symbol)+'\n')
                        log = change_log.get('rgd')
                        log[symbol] = new_val

    elif str(parser) == 'SwissWithdrawn_Parser':
        # swissprot name changes
        print('Gathering SwissProt update info...')
        cache_hits = 0
        cache_misses = 0
        files = set()
        for f in os.listdir('lions/cache/'):
            if os.path.isfile('lions/cache/'+f):
                files.add(f)
        for name in sp_lost_names:
            cached = False
            from_cache = False
            url = 'http://www.uniprot.org/uniprot/?query=mnemonic%3a'+name+ \
                '+active%3ayes&format=tab&columns=entry%20name'
            hashed_url = str(hash(url))
            ################### For Testing Only - use cache ##################
            if hashed_url in files:
                cached = True
                cache_hits += 1
                with open('lions/cache/'+hashed_url, 'rb') as fp:
                    content = pickle.load(fp)
            else:
                cache_misses += 1
                content = urllib.request.urlopen(url)

            # get the contents returned from the HTTPResponse object
            if cached:
#                ipdb.set_trace()
                content_list = content
            else:
                content_list = [x.decode().strip() for x in content.readlines()]
                with open('lions/cache/'+hashed_url, 'wb') as fp:
                    pickle.dump(content_list, fp)
            ####################################################################

            # no replacement
            if len(content_list) is 0:
                log = change_log.get('swiss-names')
                log[name] = 'withdrawn'
            # get the new name
            else:
                new_name = content_list[1]
                log = change_log.get('swiss-names')
                log[name] = new_name
        # swissprot accession numbers (ids)
        for row in parser.parse():
            sp_accession_ids.append(row.get('accession'))
        change_log['swiss-ids'] = sp_accession_ids
        print('Cache checks: ' +str(cache_hits))
        print('Cache misses: ' +str(cache_misses))

    elif str(parser) == 'GOBP_Parser':
        # GOBP name and id changes
        print('Gathering GOBP update info...')
        # treat every id as obsolete
        for lost_id in gobp_lost_ns_ids:
            log = change_log.get('gobp-ids')
            log[lost_id] = 'withdrawn'

        non_obsolete_terms = {}
        obsolete_terms = set()
        # this parser only returns the obsolete terms
        for row in parser.obsolete_parse():
            termname = row.get('termname')
            termid = row.get('termid')
            obsolete_terms.add(termid)

        # this parser only returns the NON-obsolete terms
        for row in parser.parse():
            termname = row.get('termname')
            termid = row.get('termid')
            non_obsolete_terms[termid] = termname

        # try to resolve the name
        for lost_name in gobp_lost_ns_names:
            lost_id = name_to_id(lost_name)
            if lost_id in obsolete_terms:
                log = change_log.get('gobp-names')
                log[lost_name] = 'withdrawn'
            else:
                new_name = non_obsolete_terms.get(lost_id)
                log = change_log.get('gobp-names')
                log[lost_name] = new_name

    elif str(parser) == 'GOCC_Parser':
        # GOBP name and id changes
        print('Gathering GOCC update info...')
        # treat every id as obsolete
        for lost_id in gocc_lost_ns_ids:
            log = change_log.get('gocc-ids')
            log[lost_id] = 'withdrawn'

        non_obsolete_terms = {}
        obsolete_terms = set()
        # this parser only returns the obsolete terms
        for row in parser.obsolete_parse():
            termname = row.get('termname')
            termid = row.get('termid')
            obsolete_terms.add(termid)

        # this parser only returns the NON-obsolete terms and
        # and maps them to the corresponding term-name
        for row in parser.parse():
            termname = row.get('termname')
            termid = row.get('termid')
            non_obsolete_terms[termid] = termname

        for lost_name in gocc_lost_ns_names:
            lost_id = name_to_id(lost_name)
            if lost_id in obsolete_terms:
                log = change_log.get('gocc-names')
                log[lost_name] = 'withdrawn'
            else:
                new_name = non_obsolete_terms.get(lost_id)
                log = change_log.get('gocc-names')
                log[lost_name] = new_name

    elif str(parser) == 'CHEBI_Parser':
        print('Gathering CHEBI update info...')
        # first take care of chebi, all ids are withdrawn
        log = change_log.get('chebi-ids')
        for id in chebi_lost_ids:
            log[id] = 'withdrawn'

        # some names can be resolved by using the ID
        chebi_ids_to_names = {}
        for row in parser.parse():
            new_id = row.get('primary_id')
            new_name = row.get('name')
            alt_ids = row.get('alt_ids')
            chebi_ids_to_names[new_id] = new_name
            for id in alt_ids:
                chebi_ids_to_names[id] = new_name

        log = change_log.get('chebi-names')
        # get the corresponding ID for each name
        for name in chebi_lost_names:
            lost_id = name_to_id(name)
            # if in the new namespace, use the current name as replacement
            if lost_id in new_chebi_ids:
                new_name = chebi_ids_to_names[lost_id]
                log[name] = new_name
            else:
                log[name] = 'withdrawn'

    elif str(parser) == 'MESHChanges_Parser':
        print('Gathering MESH update info...')
        old_to_new = {}
        for row in parser.parse():
            mh_old = row.get('mh_old')
            mh_new = row.get('mh_new')
            old_to_new[mh_old] = mh_new
        bio_log = change_log.get('mesh-bio')
        cell_log = change_log.get('mesh-cell')
        diseases_log = change_log.get('mesh-diseases')
        # try to resolve the lost values by using a mapping of old->new data
        for val in mesh_bio_lost:
            if val in old_to_new:
                bio_log[val] = old_to_new[val]
        for val in mesh_cell_lost:
            if val in old_to_new:
                cell_log[val] = old_to_new[val]
        for val in mesh_diseases_lost:
            if val in old_to_new:
                diseases_log[val] = old_to_new[val]

end_time = time.time()

# verification and error checking
dataset = change_log.get('entrez')
e_unresolved = [x for x in entrez_lost if x not in dataset]
dataset = change_log.get('hgnc')
hgnc_unresolved = [x for x in hgnc_lost if x not in dataset]
dataset = change_log.get('mgi')
mgi_unresolved = [x for x in mgi_lost if x not in dataset]
dataset = change_log.get('rgd')
rgd_unresolved = [x for x in rgd_lost if x not in dataset]
dataset = change_log.get('swiss-names')
sp_unresolved_names = [x for x in sp_lost_names if x not in dataset]
# values that were "lost" but not in the file provided by SwissProt
dataset = change_log.get('swiss-ids')
sp_unresolved_ids = [x for x in sp_lost_ids if x not in dataset]
dataset = change_log.get('gobp-ids')
gobp_unresolved_ids = [x for x in gobp_lost_ns_ids if x not in dataset]
dataset = change_log.get('gobp-names')
gobp_unresolved_names = [x for x in gobp_lost_ns_names if x not in dataset]
dataset = change_log.get('gocc-ids')
gocc_unresolved_ids = [x for x in gocc_lost_ns_ids if x not in dataset]
dataset = change_log.get('gocc-names')
gocc_unresolved_names = [x for x in gocc_lost_ns_names if x not in dataset]
dataset = change_log.get('chebi-names')
chebi_unresolved_names = [x for x in chebi_lost_names if x not in dataset]
dataset = change_log.get('chebi-ids')
chebi_unresolved_ids = [x for x in chebi_lost_ids if x not in dataset]
dataset = change_log.get('mesh-bio')
mesh_bio_unresolved = [x for x in mesh_bio_lost if x not in dataset]
dataset = change_log.get('mesh-cell')
mesh_cell_unresolved = [x for x in mesh_cell_lost if x not in dataset]
dataset = change_log.get('mesh-diseases')
mesh_diseases_unresolved = [x for x in mesh_diseases_lost if x not in dataset]

print('===========================================')
print('entrez unresolved: ' +str(len(e_unresolved)))
print('hgnc unresolved: ' +str(len(hgnc_unresolved)))
print('mgi unresolved: ' +str(len(mgi_unresolved)))
print('rgd unresolved: ' +str(len(rgd_unresolved)))
print('swissprot unresolved names: ' +str(len(sp_unresolved_names)))
print('swissprot unresolved ids: ' +str(len(sp_unresolved_ids)))
print('gobp unresolved names: ' +str(len(gobp_unresolved_names)))
print('gobp unresolved ids: ' +str(len(gobp_unresolved_ids)))
print('gocc unresolved names: ' +str(len(gocc_unresolved_names)))
print('gocc unresolved ids: ' +str(len(gocc_unresolved_ids)))
print('chebi unresolved names: ' +str(len(chebi_unresolved_names)))
print('chebi unresolved ids: ' +str(len(chebi_unresolved_ids)))
print('mesh-bio unresolved: ' +str(len(mesh_bio_unresolved)))
print('mesh-cell unresolved: ' +str(len(mesh_cell_unresolved)))
print('mesh-diseases unresolved: ' +str(len(mesh_diseases_unresolved)))
# add the date
today = str(datetime.date.today())
change_log['date'] = today
write_log.write(change_log)
print('total runtime is '+str(((end_time - start_time) / 60)) +' minutes.')
