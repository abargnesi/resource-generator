# coding: utf-8

'''
 datasets.py

 Represent each parsed dataset as an object. This is
 really just a wrapper to the underlying dictionaries,
 but it also provides some useful functions that assist
 in the namespacing and equivalencing process.

'''

import os.path
import time
from common import get_citation_info
from collections import defaultdict

class DataSet():
	
	def __init__(self, dictionary={}, prefix='unnamed-data-object'):
		self._dict = dictionary
		self._prefix = prefix

	def get_dictionary(self):
		return self._dict

	def __str__(self):
		return self._prefix

class OrthologyData(DataSet):

	def __init__(self, dictionary={},  prefix='use-index-term-prefix'):
		super().__init__(dictionary, prefix)
	
	def get_values(self):
		for term_id in self._dict:
			yield term_id

	def get_orthologs(self, term_id):
		orthologs = set()
		mapping = self._dict.get(term_id)
		mouse_orthologs = mapping.get('mouse_ortholog_id').split('|')
		orthologs.update(mouse_orthologs)
		if mapping.get('human_ortholog_id') is not '':
			human_orthologs = mapping.get('human_ortholog_id').split('|')
			human_orthologs = {'HGNC:' + ortho for ortho in human_orthologs}
			orthologs.update(human_orthologs)
		return orthologs

	def __str__(self):
		return self._prefix + '_ortho'

class NamespaceDataSet(DataSet):

	# Make .belns file containing ids/labels
	# default is to make .belns file for labels, and not IDs
	ids = False
	labels = True

	def __init__(self, dictionary={}, name='namespace-name', prefix='namespace-prefix'):
		self._name = name
		super().__init__(dictionary, prefix)

	def get_values(self):
		''' Get all non-obsolete primary ids in dataset dict.
		Default is all keys. '''
		for term_id in self._dict:
			yield term_id

	def get_label(self, term_id):
		''' Return the value to be used as the preferred
		label for the associated term id. Use id as default, 
		but will generally be a name/symbol. '''
		return term_id

	def get_xrefs(self, term_id):
		''' Return equivalences to other namespaces (or None). '''
		return None

	def get_name(self, term_id):
		''' Return the term name to use as title (or None). '''
		try:
			name = self._dict.get(term_id).get('name')
			return name
		except:
			return None

	def get_species(self, term_id):
		''' Return species as NCBI tax ID (or None, as applicable). '''
		return None

	def get_encoding(self, term_id):
		''' Return encoding (allowed abundance types) for value. 
		Default = 'A' (Abundance). '''
		return 'A'

	def get_alt_symbols(self, term_id):
		''' Return set of symbol synonyms. Default = None. '''
		return None

	def get_alt_names(self, term_id):
		''' Return set of name synonyms. Default = None. '''
		return None

	def get_alt_ids(self, term_id):
		''' Returns set of alternative IDs. IDs should be
		unique.  '''
		try:
			alt_ids = self._dict.get(term_id).get('alt_ids')
			return alt_ids
		except:
			return None

	def write_ns_values(self, dir):
		data_names = {}
		data_ids = {}
		for term_id in self.get_values():
			encoding = self.get_encoding(term_id)
			label = self.get_label(term_id)
			data_names[label] = encoding
			data_ids[term_id] = encoding
			if self.get_alt_ids(term_id):
				for alt_id in self.get_alt_ids(term_id):
					data_ids[alt_id] = encoding
		if self.labels:
			self.write_data(data_names, dir, self._name+'.belns')
		if self.ids:
			self.write_data(data_ids, dir, self._name+'-ids.belns')

	def write_data(self, data, dir, name):
		if len(data) == 0:
			print('	   WARNING: skipping writing ' + name + '; no namespace data found.')
		else:
			with open(os.path.join(dir, name), mode='w', encoding='utf8') as f:
				# insert header chunk
				if os.path.exists(dir+'/templates/'+name):
					tf = open(dir+'/templates/'+name, encoding="utf-8")
					header = tf.read().rstrip()
					tf.close()
					# add Namespace, Citation and Author values
					# source_file attribute added to object during parsing
					header = get_citation_info(name, header, self.source_file)
				else:
					print('WARNING: Missing header template for {0}'.format(name))
					header = '[Values]'
				f.write(header+'\n')
				# write data
				for i in sorted(data.items()):
					f.write('|'.join(i) + '\n')

	def __str__(self):
		return self._prefix

class StandardCustomData(NamespaceDataSet):
	
	def __init__(self, dictionary={}, *, name, prefix):
		super().__init__(dictionary, name, prefix)
		self._dict = {} # make unique dict for each instance of class		

	def get_values(self):
		for term_id in self._dict:
			if term_id is not None and self._dict.get(term_id).get('OBSOLETE') != 1:
				yield term_id

	def get_label(self, term_id):
		''' Return the value to be used as the preferred
		label for the associated term id. '''
		label = self._dict.get(term_id).get('LABEL')
		return label

	def get_xrefs(self, term_id):
		xrefs = self._dict.get(term_id).get('XREF').split('|')		
		return set(xrefs)
	
	def get_species(self, term_id):
		species = self._dict.get(term_id).get('SPECIES')
		return species

	def get_encoding(self, term_id):
		encoding = self._dict.get(term_id).get('TYPE')
		return encoding

	def get_alt_names(self, term_id):
		synonyms = set()
		synonyms.update(self._dict.get(term_id).get('SYNONYMS').split('|'))
		synonyms = {s for s in synonyms if s}
		return synonyms

	

class EntrezInfoData(NamespaceDataSet):

	labels = False
	ids = True
	ENC = {
		'protein-coding' : 'GRP', 'miscRNA' : 'GR', 'ncRNA' : 'GR',
		'snoRNA' : 'GR', 'snRNA' : 'GR', 'tRNA' : 'GR', 'scRNA' : 'GR',
		'other' : 'G', 'pseudo' : 'GR', 'unknown' : 'GRP', 'rRNA' : 'GR'
	}
	subject = "gene/RNA/protein"
	description = "NCBI Entrez Gene identifiers for Homo sapiens, Mus musculus, and Rattus norvegicus."

	def __init__(self, dictionary={}, name='entrez-gene', prefix='egid'):
		super().__init__(dictionary, name, prefix)

	def get_label(self, term_id):
		''' Return the value to be used as the preffered
		label for the associated term id. For Entrez, 
		using the gene ID. '''
		return term_id

	def get_species(self, term_id):
		''' Return species as NCBI tax ID (or None, as applicable). '''
		species = self._dict.get(term_id).get('tax_id')
		return species

	def get_encoding(self, gene_id):
		''' Return encoding (allowed abundance types) for value. '''
		mapping = self._dict.get(gene_id)
		gene_type = mapping.get('type_of_gene')
		description = mapping.get('description')
		encoding = EntrezInfoData.ENC.get(gene_type, 'G')
		if gene_type == 'ncRNA' and 'microRNA' in description:
			encoding = 'GRM'
		if gene_type not in EntrezInfoData.ENC:
			print('WARNING ' + gene_type + ' not defined for Entrez. G assigned as default encoding.')
		return encoding	

	def get_xrefs(self, term_id):
		''' Returns xrefs to HGNC, MGI, RGD. '''
		targets = ('MGI:', 'HGNC:', 'RGD:')
		xrefs = set()
		mapping = self._dict.get(term_id)
		xrefs.update(mapping.get('dbXrefs').split('|'))
		xrefs = {x for x in xrefs if x.startswith(targets)}
		return xrefs
			
	def get_alt_symbols(self, gene_id):
		''' Return set of symbol synonyms. '''
		synonyms = set()
		mapping = self._dict.get(gene_id)
		if mapping.get('Synonyms') is not '-':
			synonyms.update(mapping.get('Synonyms').split('|'))
			synonyms.add(mapping.get('Symbol'))
		return synonyms

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		if mapping.get('Other_designations') is not '-':
			synonyms.update(mapping.get('Other_designations').split('|'))
		if mapping.get('description') != '-':
			synonyms.add(mapping.get('description'))
		return synonyms

	def get_name(self, term_id):
		''' Get official term name. '''
		mapping = self._dict.get(term_id)
		name = mapping.get('Full_name_from_nomenclature_authority')
		return name
	

class EntrezHistoryData(DataSet):
	
	def __init__(self, dictionary={}, prefix='entrez-history'):
		super().__init__(dictionary, prefix)


class HGNCData(NamespaceDataSet, OrthologyData):
	
	ENC = {
		'gene with protein product' : 'GRP', 'RNA, cluster' : 'GR',
		'RNA, long non-coding' : 'GR', 'RNA, micro' : 'GRM',
		'RNA, ribosomal' : 'GR', 'RNA, small cytoplasmic' : 'GR',
		'RNA, small misc' : 'GR', 'RNA, small nuclear' : 'GR',
		'RNA, small nucleolar' : 'GR', 'RNA, transfer' : 'GR',
		'phenotype only' : 'G', 'RNA, pseudogene' : 'GR',
		'T cell receptor pseudogene' : 'GR',
		'immunoglobulin pseudogene' : 'GR', 'pseudogene' : 'GR',
		'T cell receptor gene' : 'GRP',
		'complex locus constituent' : 'GRP',
		'endogenous retrovirus' : 'G', 'fragile site' : 'G',
		'immunoglobulin gene' : 'GRP', 'protocadherin' : 'GRP',
		'readthrough' : 'GR', 'region' : 'G',
		'transposable element' : 'G', 'unknown' : 'GRP',
		'virus integration site' : 'G', 'RNA, micro' : 'GRM',
		'RNA, misc' : 'GR', 'RNA, Y' : 'GR', 'RNA, vault' : 'GR',
	
	}

	def __init__(self, dictionary={}, name='hgnc-human-genes', prefix='hgnc'):
		super().__init__(dictionary, name, prefix)

	def get_values(self):
		for term_id in self._dict:
			if '~withdrawn' not in self._dict.get(term_id).get('Symbol'):
				yield term_id

	def get_label(self, term_id):
		''' Return preferred label associated with term id. '''
		mapping = self._dict.get(term_id)
		if mapping is None:
			return None
		else:
			label = mapping.get('Symbol')
			return label

	def get_encoding(self, term_id):
		mapping = self._dict.get(term_id)
		locus_type = mapping.get('Locus Type')
		encoding = HGNCData.ENC.get(locus_type, 'G')
		if locus_type not in HGNCData.ENC:
			print('WARNING ' + locus_type + ' not defined for HGNC. G assigned as default encoding.')
		return encoding

	def get_species(self, term_id):
		return '9606'		

	def get_alt_symbols(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		if mapping.get('Synonyms'):
			symbol_synonyms = [s.strip() for s in mapping.get('Synonyms').split(',')]
			synonyms.update(symbol_synonyms)
		if mapping.get('Previous Symbols'):
			old_symbols = [s.strip() for s in mapping.get('Previous Symbols').split(',')]
			synonyms.update(old_symbols)
		return synonyms

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		if mapping.get('Previous Names'):
			old_names = [s.strip('" ') for s in mapping.get('Previous Names').split(', "')]
			synonyms.update(old_names)
		return synonyms

	def get_name(self, term_id):
		mapping = self._dict.get(term_id)
		name = mapping.get('Approved Name')
		return name

	def get_orthologs(self, term_id):
		orthologs = set()
		mapping = self._dict.get(term_id)
		mouse_orthologs = mapping.get('mouse_ortholog_id').split('|')
		orthologs.update(mouse_orthologs)
		rat_orthologs = mapping.get('rat_ortholog_id').split('|')
		orthologs.update(rat_orthologs)
		return orthologs

class MGIData(NamespaceDataSet):

	ENC = {
		'gene' : 'GRP', 'protein coding gene' : 'GRP',
		'non-coding RNA gene' : 'GR', 'rRNA gene' : 'GR',
		'tRNA gene' : 'GR', 'snRNA gene' : 'GR', 'snoRNA gene' : 'GR',
		'miRNA gene' : 'GRM', 'scRNA gene' : 'GR',
		'lincRNA gene' : 'GR', 'RNase P RNA gene' : 'GR',
		'RNase MRP RNA gene' : 'GR', 'telomerase RNA gene' : 'GR',
		'unclassified non-coding RNA gene' : 'GR',
		'heritable phenotypic marker' : 'G', 'gene segment' : 'G',
		'unclassified gene' : 'GRP', 'other feature types' : 'G',
		'pseudogene' : 'GR', 'transgene' : 'G',
		'other genome feature' : 'G', 'pseudogenic region' : 'GR',
		'polymorphic pseudogene' : 'GRP',
		'pseudogenic gene segment' : 'GR', 'SRP RNA gene' : 'GR'
	}

	def __init__(self, dictionary={}, name='mgi-mouse-genes', prefix='mgi'):
		super().__init__(dictionary, name, prefix)

	def get_values(self):
		for term_id in self._dict:
			mapping = self._dict.get(term_id)
			marker_type = mapping.get('Marker Type')
			if marker_type =='Gene' or marker_type == 'Pseudogene':
				yield term_id
	
	def get_species(self, term_id):
		return '10090'

	def get_encoding(self, term_id):
		feature_type = self._dict.get(term_id).get('Feature Type')
		encoding = self.ENC.get(feature_type, 'G')
		if feature_type not in self.ENC:
			print('WARNING ' + feature_type + ' not defined for MGI. G assigned as default encoding.')
		return encoding

	def get_label(self, term_id):
		try:
			label = self._dict.get(term_id).get('Symbol')
			return label
		except:
			return None

	def get_name(self, term_id):
		mapping = self._dict.get(term_id)
		name = mapping.get('Marker Name')
		return name

	def get_alt_symbols(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		synonyms = mapping.get('Marker Synonyms').split('|')
		synonyms = {s for s in synonyms if s}
		return synonyms


class RGDData(NamespaceDataSet):

	ENC = {
		'gene' : 'GRP', 'miscrna' : 'GR', 'predicted-high' : 'GRP',
		'predicted-low' : 'GRP', 'predicted-moderate' : 'GRP',
		'protein-coding' : 'GRP', 'pseudo' : 'GR', 'snrna' : 'GR',
		'trna' : 'GR', 'rrna' : 'GR', 'ncrna': 'GR'
	}

	def __init__(self, dictionary={}, name='rgd-rat-genes', prefix='rgd'):
		super().__init__(dictionary, name, prefix)
	
	def get_species(self, term_id):
		''' Rat '''
		return '10116'

	def get_label(self, term_id):
		''' Use Symbol as preferred label for RGD. '''
		try:
			label = self._dict.get(term_id).get('SYMBOL')
			return label
		except:
			return None

	def get_name(self, term_id):
		name = self._dict.get(term_id).get('NAME')
		return name

	def get_encoding(self, term_id):
		gene_type = self._dict.get(term_id).get('GENE_TYPE')
		name = self.get_name(term_id)
		encoding = RGDData.ENC.get(gene_type, 'G')
		if gene_type == 'miscrna' and 'microRNA' in name:
			encoding = 'GRM'
		if gene_type not in RGDData.ENC:
			print('WARNING ' + gene_type + ' not defined for RGD. G assigned as default encoding.')
		return encoding

	def get_alt_symbols(self, term_id):
		synonyms = set()
		if self._dict.get(term_id).get('OLD_SYMBOL'):
			old_symbols = self._dict.get(term_id).get('OLD_SYMBOL').split(';')
			synonyms.update(old_symbols)
		return synonyms

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		if mapping.get('OLD_NAME'):
			old_names = mapping.get('OLD_NAME').split(';')
			synonyms.update(old_names)
		synonyms = {s for s in synonyms if s}
		return synonyms


class SwissProtData(NamespaceDataSet):

	ids = True

	def __init__(self, dictionary=defaultdict(list), name='swissprot', prefix='sp'):
		super().__init__(dictionary, name, prefix)

	def get_encoding(self, term_id):
		return 'GRP'

	def get_label(self, term_id):
		label = self._dict.get(term_id).get('name')
		return label

	def get_name(self, term_id):
		mapping = self._dict.get(term_id)
		name = mapping.get('recommendedFullName')
		return name 

	def get_alt_ids(self, term_id):
		alt_ids = self._dict.get(term_id).get('accessions')
		alt_ids = set(alt_ids)
		alt_ids = {alt_id for alt_id in alt_ids if alt_id != term_id}
		return alt_ids

	def get_alt_symbols(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		synonyms.update(mapping.get('alternativeShortNames'))
		if mapping.get('recommendedShortName'):
			synonyms.add(mapping.get('recommendedShortname'))
		if mapping.get('geneName'):
			synonyms.add(mapping.get('geneName'))
		if mapping.get('geneSynonyms'):
			synonyms.update(mapping.get('geneSynonyms'))
		return synonyms

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		synonyms.update(mapping.get('alternativeFullNames'))
		return synonyms

	def get_xrefs(self, term_id):
		''' Returns GeneIDs or HGNC/MGI/RGD IDs. '''
		mapping = self._dict.get(term_id)
		xrefs = set()
		xrefs_dict = mapping.get('dbreference')
		for ns, values in xrefs_dict.items():
			if ns == 'GeneId':
				values = {('EGID:' + v) for v in values} 
				xrefs.update(values)
			elif ns == 'HGNC' or ns == 'MGI':
				xrefs.update(values)
			elif ns == 'RGD':
				values = {('RGD:' + v) for v in values}
				xrefs.update(values)
		return xrefs
		
	def get_species(self, term_id):
		species = self._dict.get(term_id).get('tax_id')
		return species


class AffyData(NamespaceDataSet):

	labels = False
	ids = True

	def __init__(self, dictionary=defaultdict(list), name='affy-probeset', prefix='affx'):
		super().__init__(dictionary, name, prefix)

	def get_species(self, term_id):
		species = self._dict.get(term_id).get('Species')
		species_dict = {'Homo sapiens': '9606', 
			'Mus musculus': '10090',
			'Rattus norvegicus' : '10116'}
		tax_id = species_dict.get(species)
		return tax_id

	def get_encoding(self, term_id):
		''' Return encoding (allowed abundance types) for value. 
		R - RNAAbundance. '''
		return 'R'

	def get_xrefs(self, term_id):
		''' Returns equivalent Entrez Gene IDs for value . '''
		entrez_ids = self._dict.get(term_id).get('Entrez Gene').split('///')
		if entrez_ids[0] == '---':
			return None
		else:
			entrez_ids = ['EGID:' + eid.strip() for eid in entrez_ids]
		return set(entrez_ids)	

class CHEBIData(NamespaceDataSet):

	ids = True

	def __init__(self, dictionary={}, name='chebi', prefix='chebi'):
		super().__init__(dictionary, name, prefix)
	
	def get_label(self, term_id):
		label = self._dict.get(term_id).get('name')
		return label

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		if mapping.get('synonyms'):
			synonyms.update(mapping.get('synonyms'))
		return synonyms



class Gene2AccData(DataSet):

	def __init__(self, dictionary={}, prefix = 'gene2acc'):
		super().__init__(dictionary, prefix)

	def get_eq_values(self):
		for entrez_gene in self._dict:
			mapping = self._dict.get(entrez_gene)
			status = mapping.get('status')
			taxid = mapping.get('tax_id')
			yield entrez_gene, status, taxid


class GOData(NamespaceDataSet):

	ids = True
	# dictionary is required, since GO file parsed into multiple objects
	def __init__(self, dictionary,  name, prefix):
		super().__init__(dictionary, name, prefix)

	def get_values(self):
		for term_id in self._dict:
			if self._dict.get(term_id).get('is_obsolete'):
				continue
			else:
				yield term_id

	def get_label(self, term_id):
		label = self._dict.get(term_id).get('termname')
		return label

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		synonyms.update(mapping.get('synonyms'))
		return synonyms

	def get_encoding(self, term_id):
		if self._dict.get(term_id).get('complex'):
			encoding = 'C'
		elif self._prefix == 'gobp':
			encoding = 'B'
		else:
			encoding = 'A'
		return encoding


class MESHData(NamespaceDataSet):
	# dictionary and other arguments are required since MeSH file parsed into mulitple objects
	def __init__(self, dictionary, *, name, prefix):
		super().__init__(dictionary, name, prefix)

	def get_label(self, term_id):
		label = self._dict.get(term_id).get('mesh_header')
		return label

	def get_encoding(self, term_id):
		if self._prefix == 'meshcl':
			return 'A'
		elif self._prefix == 'meshd':
			return 'O'
		elif self._prefix == 'meshpp':
			return 'B'
		else:
			return 'A'	

	def get_alt_names(self, term_id):
		synonyms = set()
		mapping = self._dict.get(term_id)
		synonyms.update(mapping.get('synonyms'))
		return synonyms	


class SwissWithdrawnData(DataSet):

	def __init__(self, dictionary={}, prefix='swiss-withdrawn'):
		super().__init__(dictionary, prefix)

	def get_withdrawn_accessions(self):
		accessions = self._dict.get('accessions')
		return accessions


class DOData(NamespaceDataSet):

	ids = True

	def __init__(self, dictionary={}, name='disease-ontology', prefix ='do'):
		super(DOData, self).__init__(dictionary, name, prefix)

	def get_label(self, term_id):
		label = self._dict.get(term_id).get('name')
		return label

	def get_encoding(self, term_id):
		return 'O'

	def get_alt_names(self,term_id):
		mapping = self._dict.get(term_id)
		synonyms  = set(mapping.get('synonyms'))
		return synonyms

	def find_xref(self, ref):
		''' Used only in equiv module. '''
		for term_id, mapping in self._dict.items():
			dbxrefs = mapping.get('dbxrefs')
			if ref in dbxrefs:
				return term_id

	def get_xrefs(self, term_id):
		''' Returns MeSH (MSH) xrefs for a given DO ID . '''
		xrefs = set()
		mapping = self._dict.get(term_id)
		xrefs.update(mapping.get('dbxrefs'))
		xrefs = {x.replace('MSH:','MESHD:') for x in xrefs if x.startswith('MSH:')}
		return xrefs
# vim: ts=4 sts=4 sw=4 noexpandtab
