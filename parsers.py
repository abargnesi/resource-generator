# coding: utf-8

'''
 parsers.py

 Contains parser objects constructed for various datasets that
 will be used to build the .belns, .beleq, and .belanno files.

'''

from common import gzip_to_text
from lxml import etree
from collections import defaultdict
import os
import csv
import gzip
import urllib.request
import zipfile
import io

class Parser(object):
	''' Generic/parent parser. '''
		
	def __init__(self, url):
		self._url = url
		self.verbose = False

	def is_verbose(self):
		self.verbose = True

	def parse(self):
		with open(self._url) as f:
			reader = csv.DictReader(f, delimiter='\t')
			for row in reader:
				yield row
	
	def __str__(self):
		return "Parser"

class RGDOrthologParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		with open(self._url, 'r') as f:
			reader = csv.DictReader(filter(lambda row: 
								not row[0].startswith('#'), f), delimiter='\t')
			for row in reader:
				yield row

	def __str__(self):
		return "RGDOrthologParser"

class NamespaceParser(Parser):
	''' Generic parser. Expects tab-delimited file - 
	all rows prior to column header row should start with '#'.  '''
		
	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		with open(self._url, 'r') as f:
			reader = csv.DictReader(filter(lambda row: 
								not row[0].startswith('#'), f), delimiter='\t')
			for row in reader:
				if row['ID'].strip():
					yield row

	def __str__(self):
		return "NamespaceParser"

class EntrezGeneInfoParser(Parser):

	headers = ['tax_id', 'GeneID', 'Symbol', 'LocusTag',
					   'Synonyms', 'dbXrefs', 'chromosome',
					   'map_location', 'description',
					   'type_of_gene',
					   'Symbol_from_nomenclature_authority',
					   'Full_name_from_nomenclature_authority',
					   'Nomenclature_status',
					   'Other_designations', 'Modification_date']

	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		reader = csv.DictReader(gzip_to_text(self._url),
					   delimiter='\t',
					   fieldnames=self.headers)

		for row in reader:
			if row['tax_id'] in ('9606', '10090', '10116'):
				yield row
            #    continue     

	def __str__(self):
		return "EntrezGeneInfo_Parser"


class EntrezGeneHistoryParser(Parser):

	headers = ["tax_id", "GeneID", "Discontinued_GeneID",
							  "Discontinued_Symbol", "Discontinue_Date"]

	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		reader = csv.DictReader(gzip_to_text(self._url),
						  delimiter='\t',
						  fieldnames=self.headers)

		for row in reader:
			if row['tax_id'] in ("9606", "10090", "10116"):
				yield row

	def __str__(self):
		return "EntrezGeneHistory_Parser"


class HGNCParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):

		# use iso-8859-1 as default encoding.
		with open(self._url, "r", encoding="iso-8859-1") as hgncf:

			# Note that HGNC uses TWO columns named the same thing for Entrez
			# Gene ID. Currently we are not using these columns and it is not a
			# big deal, but in the future we could account for this by using
			# custom headers (like EntrezGeneInfo_Parser), or resolving to the
			# SECOND of the Entrez Gene ID columns.
			reader = csv.DictReader(hgncf, delimiter='\t')

			for row in reader:
				yield row

	def __str__(self):
		return "HGNC_Parser"


class MGIParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		with open(self._url, "r") as f:
			reader = csv.DictReader(f, delimiter='\t')

			for row in reader:
				yield row

	def __str__(self):
		return "MGI_Parser"


class RGDParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):
		with open(self._url, "r") as f:
			# skip all the comment lines beginning with '#' and also the header.
			reader = csv.DictReader(filter(lambda row:
									 not row[0].startswith('#'), f),
									  delimiter='\t')

			for row in reader:
				yield row

	def __str__(self):
		return "RGD_Parser"


class GeneTypeError(Exception):
	''' This class exists mainly as a way to break the iteration loop during parsing
	 of the SwissProt dataset if needed.'''
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class SwissProtParser(Parser):

	def __init__(self, url):
		super().__init__(url)
		self.entries = {}
		self.accession_numbers = {}
		self.gene_ids = {}
		self.tax_ids = {'9606', '10090', '10116'}
		self.pro = '{http://uniprot.org/uniprot}protein'
		self.rec_name = '{http://uniprot.org/uniprot}recommendedName'
		self.full_name = '{http://uniprot.org/uniprot}fullName'
		self.short_name = '{http://uniprot.org/uniprot}shortName'
		self.alt_name = '{http://uniprot.org/uniprot}alternativeName'
		self.db_ref = '{http://uniprot.org/uniprot}dbReference'
		self.organism = '{http://uniprot.org/uniprot}organism'
		self.entry = '{http://uniprot.org/uniprot}entry'
		self.accession = '{http://uniprot.org/uniprot}accession'
		self.name = '{http://uniprot.org/uniprot}name'
		self.gene = '{http://uniprot.org/uniprot}gene'

	def parse(self):

		with gzip.open(self._url) as f:
			ctx = etree.iterparse(f, tag=self.entry)

			for ev, e in ctx:
				temp_dict = {}
				n_dict = defaultdict(list)

				# stop evaluating if this entry is not in the Swiss-Prot dataset
				if e.get('dataset') != 'Swiss-Prot':
					e.clear()
					continue

				# stop evaluating if this entry is not for human, mouse, or rat
				org = e.find(self.organism)

				# use a custom exception to break to next iteration (e)
				# if tax ref is not found.
				try:
					for org_child in org:
						if org_child.tag == self.db_ref:
							# restrict by NCBI Taxonomy reference
							if org_child.get('id') not in self.tax_ids:
								e.clear()
								raise GeneTypeError(org_child.get('id'))
							else:
								# add NCBI Taxonomy and the id for the entry
								# to the dict
								temp_dict[org_child.get('type')] = \
									org_child.get('id')
				except GeneTypeError:
					continue

				# get entry name, add it to the dict
				entry_name = e.find(self.name).text
				temp_dict['name'] = entry_name
 
				# get protein data, add recommended full and short names to dict
				protein = e.find(self.pro)

				for child in protein.find(self.rec_name):
					if child.tag == self.full_name:
						temp_dict['recommendedFullName'] = child.text
					if child.tag == self.short_name:
						temp_dict['recommendedShortName'] = child.text
				alt_shortnames = []
				alt_fullnames = []

				protein = e.find(self.pro)
				for altName in protein.findall(self.alt_name):
					for child in altName:
						if child.tag == self.full_name:
							alt_fullnames.append(child.text)
						if child.tag == self.short_name:
							alt_shortnames.append(child.text)

				temp_dict['alternativeFullNames'] = alt_fullnames
				temp_dict['alternativeShortNames'] = alt_shortnames

				# get gene data, add primary names (symbols) and synonyms
				gene = e.find(self.gene)
				gene_synonyms = []
				if gene is not None:
					for name in gene.findall(self.name):
						if name.get('type') == 'primary':
							gene_name = name.text
						elif name.get('type') == 'synonym':
							gene_synonyms.append(name.text)
					temp_dict['geneName'] = gene_name
					temp_dict['geneSynonyms'] = gene_synonyms
				
				# get all accessions
				entry_accessions = []
				for entry_accession in e.findall(self.accession):
					acc = entry_accession.text
					entry_accessions.append(acc)
					if acc in self.accession_numbers:
						self.accession_numbers[acc] = None
					else:
						self.accession_numbers[acc] = 1

				# add the array of accessions to the dict
				temp_dict["accessions"] = entry_accessions

				# add dbReference type (human, rat, and mouse) and gene ids to
				# the dict
				type_set = ['GeneId', 'MGI', 'HGNC', 'RGD']
				for dbr in e.findall(self.db_ref):
					if dbr.get('type') in type_set:
						gene_id = dbr.get('id')
						n_dict[dbr.get('type')].append(gene_id)
				temp_dict['dbreference'] = n_dict

				# clear the tree before next iteration
				e.clear()
				while e.getprevious() is not None:
					del e.getparent()[0]
				yield temp_dict

	def __str__(self):
		return 'SwissProt_Parser'


# Helper function for AffyParser. This will save each of the downloaded
# URLs and return the file pointer.
def get_data(url):
	# from url, download and save file
	REQ = urllib.request.urlopen(url)
	file_name = url.split('/')[-1]
	os.chdir('datasets/')
	with open(file_name,'wb') as f:
		f.write(REQ.read())
	os.chdir('../')
	return file_name

def filter_plus_print(row):
	return not row.startswith('#')


class AffyParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):

		from configuration import affy_array_names
		urls = []
		with open(self._url, 'rb') as f:
			ctx = etree.iterparse(f, events=('start', 'end'))

			# This is probably not the best way to traverse this tree. Look at
			# the lxml.etree API more closely for possible implementations when
			# refactoring
			# NOTES - put some debugging in here to see how this is parsing,
			# may be a better way to parse (like using diff events).
			for ev, e in ctx:
				# iterate the Array elements
				for n in e.findall('Array'):
					name = n.get('name')
					if name in affy_array_names:
						# iterate Annotation elements
						for child in n:
							if child.get('type') == 'Annot CSV':
								# iterate File elements
								for g_child in child:
									# get the URL and add to the list
									for gg_child in g_child:
										urls.append(gg_child.text)

		# iterate over the list of URLs returned from the Affy XML feed
		for link in urls:
			affy_reader = {}

			# get_data() downloads the file, saves it as a .csv.zip, and
			# returns a pointer to the file.
			n = get_data(link)
			z = zipfile.ZipFile('datasets/'+n, 'r')

			# only want the .csv from the archive (also contains a .txt)
			for name in z.namelist():
				if '.csv' in name:
					if self.verbose:
						print('\tExtracting - ' +name)
					# wrap in a TextIOWrapper. otherwise it returns bytes.
					affy_reader = csv.DictReader(filter(lambda x:
														not x.startswith('#'),
														io.TextIOWrapper(z.open(name))),
														delimiter=',')

					for x in affy_reader:
						yield x

	def __str__(self):
		return 'Affy_Parser'

class Gene2AccParser(Parser):

	def __init__(self, url):
	   super().__init__(url)

	def parse(self):

		# would like to have DictReader handle this, but need a way to
		# deal with the special case of the first value beginning with
		# a hashtag. i.e. #Format: <-- is NOT a column header.
		column_headers = ['tax_id', 'GeneID', 'status',
						  'RNA nucleotide accession.version',
						  'RNA nucleotide gi', 'protein accession.version',
						  'protein gi', 'genomic nucleotide accession.version',
						  'genomic nucleotide gi',
						  'start position on the genomic accession',
						  'end position on the genomic accession',
						  'orientation', 'assembly',
						  'mature peptide accession.version',
						  'mature peptide gi', 'Symbol']

		g2a_reader = csv.DictReader(gzip_to_text(self._url), delimiter='\t',
									fieldnames=column_headers)

		for row in g2a_reader:
			if row['tax_id'] in ('9606', '10090', '10116'):
				yield row

	def __str__(self):
		return 'Gene2Acc_Parser'

class BELNamespaceParser(Parser):

	def __init__(self):
		self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'
		self.anno_def = '{http://www.belscript.org/schema/annotationdefinitions}annotationdefinitions'
		self.namespace = '{http://www.belscript.org/schema/namespace}namespace'
		self.namespaces = '{http://www.belscript.org/schema/namespaces}namespaces'

	def parse(self):

		tree = etree.parse(self.old_files)

		# xpath will return all elements under this namespace (list of bel namespace urls)
		urls = tree.xpath('//*[local-name()="namespace"]/@idx:resourceLocation',
						  namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
									  'idx' : 'http://www.belscript.org/schema/index'})

		for url in urls:
			yield url

	def __str__(self):
		return 'BELNamespace_Parser'


class BELEquivalenceParser(Parser):

	def __init__(self):
		self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'
		self.anno_def = '{http://www.belscript.org/schema/annotationdefinitions}annotationdefinitions'
		self.namespace = '{http://www.belscript.org/schema/namespace}namespace'
		self.namespaces = '{http://www.belscript.org/schema/namespaces}namespaces'

	def parse(self):

		tree = etree.parse(self.old_files)

		# xpath will return all elements under this namespace (list of bel equivalence urls)
		urls = tree.xpath('//*[local-name()="equivalence"]/@idx:resourceLocation',
						  namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
									  'idx' : 'http://www.belscript.org/schema/index'})

		for url in urls:
			yield url

	def __str__(self):
		return 'BELEquivalence_Parser'


class BELAnnotationsParser(Parser):

	def __init__(self):
		self.old_files = 'http://resource.belframework.org./belframework/1.0/index.xml'

	def parse(self):

		tree = etree.parse(self.old_files)

		# xpath will return all elements under this namespace (list of bel equivalence urls)
		urls = tree.xpath('//*[local-name()="annotationdefinition"]/@idx:resourceLocation',
						  namespaces={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
									  'idx' : 'http://www.belscript.org/schema/index'})

		for url in urls:
			yield url

	def __str__(self):
		return 'BELAnnotations_Parser'

# This one uses iterparse(), much faster than xpath on the
# bigger .owl file.
class CHEBIParser(Parser):

	def __init__(self, url):
		super().__init__(url)
		self.classy = '{http://www.w3.org/2002/07/owl#}Class'
		self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'
		self.altId = '{http://purl.obolibrary.org/obo#}altId'
		self.synonym = '{http://purl.obolibrary.org/obo#}Synonym'

	def parse(self):

		with open(self._url, 'rb') as cf:
			tree = etree.iterparse(cf, tag=self.classy)
			for event, elem in tree:
				if len(elem.values()) != 0:
					chebi_dict = {}
					synonyms = set()
					alt_ids = set()
					name = ''
					vals = elem.values()
					chebi_dict['primary_id'] = vals[0].split('CHEBI_')[1]
					children = elem.getchildren()
					for child in children:
						if child.tag == self.label:
							name = child.text
						if child.tag == self.altId:
							alt_ids.add(child.text.split(':')[1])
						if child.tag == self.synonym:
							synonyms.add(child.text)
						chebi_dict['name'] = name
						chebi_dict['alt_ids'] = alt_ids
						chebi_dict['synonyms'] = synonyms

					yield chebi_dict

	def __str__(self):
		return 'CHEBI_Parser'



class GOParser(Parser):
	# TODO - build obsolete methods into data set
	def __init__(self, url):
		super().__init__(url)

	def parse(self):

		# initialize empty dictionaries using tuple assignment
		parents, accession_dict, term_dict = {}, {}, {}

		# parse xml tree using lxml
		parser = etree.XMLParser(ns_clean=True, recover=True, encoding='UTF-8')
		root = etree.parse(self._url, parser)
		terms = root.xpath("/obo/term")

		# iterate the complex terms to build parent dictionary
		for t in terms:
			if t.find('namespace').text == 'cellular_component':
				termid = t.find("id").text
				parent_ids = [isa.text for isa in t.findall("is_a")]
				parents[termid] = parent_ids

		for t in terms:
			is_obsolete = False
			termid = t.find('id').text
			termname = t.find('name').text
			namespace = t.find('namespace').text
			if t.find('is_obsolete') is not None:
				is_obsolete = True
			if t.findall('alt_id') is not None:
				altids = [x.text for x in t.findall('alt_id')]
			else:
				altids = False
	
			# identify complexes (for GOCC), based on parent terms
			is_complex = False
			if namespace == 'cellular_component':
				parent_stack = parents[termid]
				if termid == "GO:0032991":
					is_complex = True
				elif t.find("is_root") is not None:
					is_complex = False
				else:
					parent_stack.extend(parents[termid])
					while len(parent_stack) > 0:
						parent_id = parent_stack.pop()
						if parent_id == "GO:0032991":
							is_complex = True
							break
						if parent_id in parents:
							parent_stack.extend(parents[parent_id])

			# strip 'GO:' from term_ids	
			termid = termid.replace('GO:','')
			altids = [altid.replace('GO:','') for altid in altids]
 
			#get synonyms - limited to scope='exact'
			synonyms = []
			for syn in t.findall('synonym'):
				if syn.get('scope') == 'exact':
					synonyms.append(syn.find('synonym_text').text)

			yield { 'termid' : termid, 'termname' : termname, 'namespace' : namespace,
					'altids' : altids, 'complex' : is_complex, 'synonyms' : synonyms,
					'is_obsolete' : is_obsolete }

	def __str__(self):
		return 'GO_Parser'

class MESHParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):

		# ui - unique identifier / mh - mesh header
		# mn - tree # / st - semantic type
		ui = ''
		mh = ''
		mns = set()
		sts = set()
		synonyms = set()
		firstTime = True
		with open(self._url, 'r') as fp:
			for line in fp.readlines():
				if line.startswith('MH ='):
					mh = line.lstrip('MH =').strip()
				elif line.startswith('UI ='):
					ui = line.lstrip('UI =').strip()
				elif line.startswith('MN ='):
					mn = line.lstrip('MN =').strip()
					mns.add(mn)
				elif line.startswith('ST ='):
					st = line.lstrip('ST =').strip()
					sts.add(st)
				elif line.startswith('PRINT ENTRY ='):
					entry = line.lstrip('PRINT ENTRY =').strip()
					if '|EQV|' in line:
						entries = entry.split('|')
						last = entries[-1]
						num_syns = last.count('a')
						while num_syns > 0:
							num_syns = num_syns - 1
							s = entries[num_syns]
							synonyms.add(s.strip())
					else:
						if '|' not in entry:
							synonyms.add(entry)
				elif line.startswith('ENTRY ='):
					entry = line.lstrip('ENTRY =').strip()
					if '|EQV|' in line:
						entries = entry.split('|')
						last = entries[-1]
						num_syns = last.count('a')
						while num_syns > 0:
							num_syns = num_syns - 1
							s = entries[num_syns]
							synonyms.add(s.strip())
					else:
						if '|' not in entry:
							synonyms.add(entry)
				elif line.startswith('*NEWRECORD'):
					# file begins with *NEWRECORD so skip that one (dont yield)
					if firstTime:
						firstTime = False
						continue
					else:

						yield { 'ui' : ui, 'mesh_header' : mh,
								'mns' : mns, 'sts' : sts,
								'synonyms' : synonyms }
						ui = ''
						mh = ''
						mns = set()
						sts = set()
						synonyms = set()

	def __str__(self):
		return 'MESH_Parser'


class SwissWithdrawnParser(Parser):

	def __init__(self, url):
		super(SwissWithdrawnParser, self).__init__(url)
		self.s_file = url

	def parse(self):

		with open(self.s_file, 'r') as fp:
			marker = False
			for line in fp.readlines():
				if '____' in line:
					marker = True
					continue
				if marker is False:
					continue

				yield {'accession' : line.strip()}

	def __str__(self):
		return 'SwissWithdrawn_Parser'


class MESHChangesParser(Parser):

	def __init__(self, url):
		super(MESHChangesParser, self).__init__(url)
		self.mesh_file = url

	def parse(self):

		with open(self.mesh_file, 'r') as fp:
			for line in fp.readlines():
				if 'MH OLD =' in line:
					mh_old = line.split('= ')[1]
					if '#' in mh_old:
						mh_old = mh_old.split(' #')[0]
					elif '[' in mh_old:
						mh_old = mh_old.split(' [')[0]
				if 'MH NEW =' in line:
					mh_new = line.split('= ')[1]
					if '#' in mh_new:
						mh_new = mh_new.split(' #')[0]
					elif '[' in mh_new:
						mh_new = mh_new.split(' [')[0]
					yield { 'mh_old' : mh_old.strip(), 'mh_new' : mh_new.strip() }
					mh_old = ''
					mh_new = ''

	def __str__(self):
		return 'MESHChanges_Parser'


def is_deprecated(child_list):
	dep = False
	deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
	for child in child_list:
		if child.tag == deprecated and child.text == 'true':
			dep = True
	return dep


# custom exception to break out of the loop when an deprecated
# term has been seen.
class DeprecatedTermException(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

# excludes deprecated terms which are not included in the namespace
class DOParser(Parser):

	def __init__(self, url):
		super().__init__(url)
		self.classy = '{http://www.w3.org/2002/07/owl#}Class'
		self.deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
		self.dbxref = '{http://www.geneontology.org/formats/oboInOwl#}hasDbXref'
		self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'
		self.exactsynonym = '{http://www.geneontology.org/formats/oboInOwl#}hasExactSynonym'
		about = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'

	def parse(self):

		with open(self._url, 'rb') as df:
			tree = etree.iterparse(df, tag=self.classy)
			about = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about'
			for event, elem in tree:
				do_dict = {}
				dbxrefs = []
				synonyms = []
				name = ''
				term_id = ''
				try:
					if len(elem.values()) != 0:
						children = elem.getchildren()
						if is_deprecated(children):
							raise DeprecatedTermException(children)
						else:
							term_id = elem.get(about).split('/')[-1].strip('DOID_')
							for child in children:
								if child.tag == self.dbxref:
									dbxrefs.append(child.text)
								elif child.tag == self.label:
									name = child.text
								elif child.tag == self.exactsynonym:
									synonyms.append(child.text)
							do_dict['name'] = name
							do_dict['id'] = term_id
							do_dict['dbxrefs'] = dbxrefs
							do_dict['synonyms'] = synonyms
							yield do_dict
				except DeprecatedTermException:
					continue

	def __str__(self):
		return 'DO_Parser'


# includes deprecated terms (for change-log)
class DODeprecatedParser(Parser):

	def __init__(self, url):
		super().__init__(url)
		self.classy = '{http://www.w3.org/2002/07/owl#}Class'
		self.deprecated = '{http://www.w3.org/2002/07/owl#}deprecated'
		self.dbxref = '{http://www.geneontology.org/formats/oboInOwl#}hasDbXref'
		self.id = '{http://www.geneontology.org/formats/oboInOwl#}id'
		self.label = '{http://www.w3.org/2000/01/rdf-schema#}label'

	def parse(self):

		with open(self._url, 'rb') as df:
			tree = etree.iterparse(df, tag=self.classy)
			for event, elem in tree:
				do_dep_dict = {}
				dbxrefs = []
				name = ''
				id = ''
				dep = False
				if len(elem.values()) != 0:
					children = elem.getchildren()
					for child in children:
						if child.tag == self.dbxref:
							dbxrefs.append(child.text)
						elif child.tag == self.id:
							id = child.text.split(':')[1]
						elif child.tag == self.label:
							name = child.text
						elif child.tag == self.deprecated:
							dep = True
				do_dep_dict['name'] = name
				do_dep_dict['id'] = id
				do_dep_dict['dbxrefs'] = dbxrefs
				do_dep_dict['deprecated'] = dep
				yield do_dep_dict

	def __str__(self):
		return 'DODeprecated_Parser'

#TODO - consolidate RGD parsers - file format is same
# filtering of non-rat data can be done in parsed module
class RGDObsoleteParser(Parser):

	def __init__(self, url):
		super().__init__(url)

	def parse(self):

		with open(self._url, 'r') as rgdo:
			# skip comment lines
			rgd_csvr = csv.DictReader(filter(lambda row:
												not row[0].startswith('#'), rgdo), 
									 delimiter='\t')
			for row in rgd_csvr:
				if row['SPECIES'] == 'rat':
					yield row

	def __str__(self):
		return "RGD_Obsolete_Parser"

# vim: ts=4 sts=4 sw=4 noexpandtab
