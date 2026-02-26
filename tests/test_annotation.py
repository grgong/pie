from pie.annotation import parse_annotations, GeneModel


class TestParseAnnotations:
    def test_parse_gff3(self, gff3_file):
        genes = parse_annotations(gff3_file)
        assert len(genes) == 3

    def test_parse_gtf(self, gtf_file):
        genes = parse_annotations(gtf_file)
        assert len(genes) == 3

    def test_same_gene_count(self, gff3_file, gtf_file):
        assert len(parse_annotations(gff3_file)) == len(parse_annotations(gtf_file))


class TestGeneModel:
    def test_gene1_attributes(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g1 = genes[0]
        assert "gene1" in g1.gene_id.lower()
        assert g1.strand == "+"
        assert g1.cds_length == 90
        assert len(g1.cds_exons) == 1

    def test_gene2_multi_exon(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g2 = [g for g in genes if "gene2" in g.gene_id.lower()][0]
        assert len(g2.cds_exons) == 2
        assert g2.cds_length == 99  # 59 + 40
        assert g2.strand == "+"

    def test_gene3_minus_strand(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g3 = [g for g in genes if "gene3" in g.gene_id.lower()][0]
        assert g3.strand == "-"
        assert g3.cds_length == 81
        assert len(g3.cds_exons) == 1

    def test_cds_exons_sorted(self, gff3_file):
        genes = parse_annotations(gff3_file)
        for gene in genes:
            starts = [s for _, s, _ in gene.cds_exons]
            assert starts == sorted(starts)

    def test_coordinates_0based(self, gff3_file):
        genes = parse_annotations(gff3_file)
        g1 = genes[0]
        # GFF says 1-90, should be 0-based: start=0, end=90 (half-open)
        assert g1.cds_exons[0] == ("chr1", 0, 90)
