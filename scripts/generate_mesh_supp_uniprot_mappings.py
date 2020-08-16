"""Generate MeSH supplementary concept -> Uniprot Mappings

This script identifies MeSH supplementary concept entries corrresponding to human
proteins. It works by finding supplementary concepts with name or a synonym of the
form 'X protein, human' where X is a valid HGNC symbol. If a supplementary concept
has a unique such X then the match is considered to be of type 'skos:exactMatch.
If there are multiple such Xs then the match is considered to be of type
'skos:narrowMatch'.
"""


import re
import os
import gzip
import urllib
from lxml import etree
from subprocess import check_output
from indra.databases import mesh_client, hgnc_client


MESH_SUPP_URL = ('ftp://nlmpubs.nlm.nih.gov/'
                 'online/mesh/MESH_FILES/xmlmesh/supp2020.gz')

def get_script_url():
    commit_hash = check_output('git rev-parse HEAD'.split()).\
        decode('utf-8').strip()[:6]
    script_name = os.path.basename(__file__)
    return (f'https://github.com/biomappings/biomappings/blob/{commit_hash}/'
            f'scripts/{script_name}')


def get_mappings():
    url = get_script_url()
    mapping_type = 'lexical'
    ftp_request = urllib.request.Request(MESH_SUPP_URL)
    with urllib.request.urlopen(ftp_request) as fp:
        with gzip.GzipFile(fileobj=fp) as gz_file:
            tree = etree.parse(gz_file)
    for concept_element in tree.xpath('SupplementalRecord'):
        UI_element = concept_element.xpath('SupplementalRecordUI')
        if not UI_element:
            continue
        mesh_id = UI_element[0].text
        name_element = concept_element.xpath('SupplementalRecordName/String')
        if not name_element:
            continue
        mesh_name = name_element[0].text
        synonyms = [element.text for
                    element in concept_element.xpath('.//TermList/Term/String')]
        matched_genes = set()
        for term in [mesh_name] + synonyms:
            match = re.match(r'^(.+) protein, human$', term)
            if match:
                gene_name = match.groups()[0]
                hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                if hgnc_id is not None:
                    matched_genes.add((gene_name, hgnc_id))
        if len(matched_genes) == 1:
            match_type = 'skos:exactMatch'
        else:
            match_type = 'skos:narrowMatch'
        for gene_name, hgnc_id in matched_genes:
            uniprot_id = hgnc_client.get_uniprot_id(hgnc_id)
            if uniprot_id is not None:
                yield ('mesh', mesh_id, mesh_name, match_type, 'uniprot',
                       uniprot_id, gene_name, mapping_type, url)
