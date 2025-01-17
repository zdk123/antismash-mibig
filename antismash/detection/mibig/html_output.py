# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Handles HTML output for the MIBiG sideloader """

import os
from typing import List

from eutils import Client

from antismash.common import path
from antismash.common.module_results import ModuleResults
from antismash.common.html_renderer import FileTemplate, HTMLSections
from antismash.common.layers import RegionLayer, RecordLayer, OptionsLayer

from mibig.converters.read.cluster import Publication

from .mibig import MibigAnnotations


def will_handle(_products: List[str]) -> bool:
    """ Returns true if one or more relevant products are present """
    return True


def generate_html(region_layer: RegionLayer, results: ModuleResults,
                  record_layer: RecordLayer, _options_layer: OptionsLayer) -> HTMLSections:
    assert isinstance(results, MibigAnnotations)
    data = results.data
    taxonomy = results.taxonomy

    html = HTMLSections("mibig-general")
    taxonomy_text = " > ".join(["{} ({})".format(taxobj["name"], taxobj["rank"]) for taxobj in taxonomy])
    publications_links = ReferenceCollection(data.cluster.publications)
    html.add_detail_section("General", FileTemplate(path.get_full_path(__file__, "templates", "general.html")).render(data=results.data, taxonomy_text=taxonomy_text, publications_links=publications_links.links()))

    html.add_detail_section("Compounds", FileTemplate(path.get_full_path(__file__, "templates", "compounds.html")).render(compounds=results.data.cluster.compounds),
                            class_name="mibig-compounds")

    genes = []
    annots = data.cluster.genes.annotations if data.cluster.genes else []
    for cds_feature in record_layer.get_cds_features_within_location(region_layer.location):
        gene = {
            "locus_tag": cds_feature.locus_tag,
            "protein_id": cds_feature.protein_id,
            "gene": cds_feature.gene,
            "start": cds_feature.location.start + 1,
            "end": cds_feature.location.end,
            "strand": "+" if cds_feature.location.strand > 0 else "-",
            "product": cds_feature.product,
            "aa_seq": cds_feature.translation,
            "nt_seq": cds_feature.extract(record_layer.seq)
        }
        gene["functions"] = []
        for function in cds_feature.gene_functions:
            function_text = str(function.function)
            if function.tool != "mibig":
                continue
            if function.tool == "rule-based-clusters":
                function_text += " ({})".format(function.description)
            elif function.tool == "smcogs":
                function_text += " ({})".format(function.description.split(" (")[0])
            gene["functions"].append(function_text)
        annot_idx = -1
        for i, annot in enumerate(annots):
            name = annot.name
            if annot.id == cds_feature.locus_tag or annot.id == cds_feature.protein_id or (name and name == cds_feature.gene):
                annot_idx = i
                break
        if annot_idx >= 0:
            annot = annots.pop(annot_idx)
            for function in annot.functions:
                function_text = function.category
                if annot.tailoring:
                    function_text += " ({}) ".format(", ".join(annot.tailoring))
                else:
                    function_text += " "
                gene["functions"].append(function_text)
                gene["evidences"] = sorted(set(function.evidence))
            if annot.mutation_phenotype:
                function_text = "Mutation phenotype: {}".format(annot.mutation_phenotype)
                gene["functions"].append(function_text)
            if annot.product:
                gene["product"] = annot.product
        genes.append(gene)
    for annot in annots:
        gene = {
            "locus_tag": annot.id,
            "protein_id": "None",
            "gene": annot.name,
            "product": annot.product or ""
        }
        gene["functions"] = []
        for function in annot.functions:
            function_text = function.category
            if annot.tailoring:
                function_text += " ({}) ".format(", ".join(annot.tailoring))
            else:
                function_text += " "
            gene["functions"].append(function_text)
            gene["evidences"] = sorted(set(function.evidence))
        genes.append(gene)
    html.add_detail_section("Genes", FileTemplate(path.get_full_path(__file__, "templates", "genes.html")).render(genes=genes),
                            class_name="mibig-genes")

    if data.cluster.polyketide:
        html.add_detail_section("Polyketide", FileTemplate(path.get_full_path(__file__, "polyketide", "details.html")).render(pk=results.data.cluster.polyketide),
                                class_name="mibig-polyketide")

    if data.cluster.nrp:
        html.add_detail_section("NRP", FileTemplate(path.get_full_path(__file__, "nrp", "details.html")).render(nrp=results.data.cluster.nrp),
                                class_name="mibig-nrp")

    if data.cluster.ripp:
        html.add_detail_section("RiPP", FileTemplate(path.get_full_path(__file__, "ripp", "details.html")).render(ripp=results.data.cluster.ripp),
                                class_name="mibig-ripp")

    if data.cluster.saccharide:
        html.add_detail_section("Saccharide", FileTemplate(path.get_full_path(__file__, "saccharide", "details.html")).render(sac=results.data.cluster.saccharide),
                                class_name="mibig-saccharide")

    if data.cluster.terpene:
        html.add_detail_section("Terpene", FileTemplate(path.get_full_path(__file__, "terpene", "details.html")).render(trp=results.data.cluster.terpene),
                                class_name="mibig-terpene")

    html.add_detail_section("History", FileTemplate(path.get_full_path(__file__, "templates", "logs.html")).render(logs=sorted(results.data.changelog, key=lambda log: log.mibig_version)),
                                class_name="mibig-logs")

    return html


class ReferenceLink:
    """Keep track of a single reference link."""

    __slots__ = (
        'category',
        'ref',
        'title',
        'info',
    )

    def __init__(self, category: str, reference: str, title: str, info: str = None) -> None:
        self.category = category
        self.ref = reference
        self.title = title
        self.info = info


class ReferenceCollection:
    """Keep track of all references in a MIBiG entry."""

    __slots__ = (
        'client',
        'references',
    )

    def __init__(self, publications: List[Publication]) -> None:
        self.client = Client(api_key=os.environ.get("NCBI_API_KEY", None))
        self.references = {}  # Dict[ReferenceLink]
        pmids = []

        for publication in publications:
            info = None

            if publication.category == "pubmed":
                reference = "https://www.ncbi.nlm.nih.gov/pubmed/{}".format(publication.content)
                # TODO: Remove once BGC0001346 is fixed
                if publication.content != "0":
                    pmids.append(publication.content)
            elif publication.category == "patent":
                reference = "https://patents.google.com/patent/{}".format(publication.content)
            elif publication.category == "doi":
                reference = "https://dx.doi.org/{}".format(publication.content)
            elif publication.category == "url":
                reference = publication.content

            self.references[publication.content] = ReferenceLink(publication.category, reference, publication.content, info)

        self._resolve_pmids(pmids)


    def links(self) -> List[ReferenceLink]:
        return self.references.values()


    def _resolve_pmids(self, pmids: List[str]) -> None:
        if not pmids:
            return
        articles = self.client.efetch(db="pubmed", id=pmids)
        for article in articles:
            self.references[article.pmid].title = article.title
            self.references[article.pmid].info = "{a.authors[0]} et al., {a.jrnl} ({a.year}) PMID:{a.pmid}".format(a=article)
